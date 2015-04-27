
using ArgParse

function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--input", "-i"
      help = "input FASTQ filename"
    "--barcode"
      help = "barcode file, tab delimited"
      arg_type = ASCIIString
    "--prefix"
      help = "prefix for split files"
      arg_type = ASCIIString
      default = ""
    "--suffix"
      help = "suffix for split files"
      arg_type = ASCIIString
      default = ".fastq"
    "--keep_unmatched"
      help = "save reads unmatched to any barcodes"
      action = :store_true
    "--mismatches", "-m"
      help = "maximum number of mismatches"
      arg_type = Int
      default = 0
    "--partial", "-p"
      help = "number of partial matches allowed"
      arg_type = Int
      default = 0
  end

  return parse_args(s)
end

parsed_args = parse_commandline()

input = open(parsed_args["input"],"r")

MISMATCHES = parsed_args["mismatches"]
PARTIAL = parsed_args["partial"]
KEEP_UNMATCHED = parsed_args["keep_unmatched"]

code_map = convert(Array{ASCIIString},readdlm(parsed_args["barcode"]))
n_codes = size(code_map)[1]

code_files = Array(IOStream,n_codes+1)
code_counts = Array(Int64,n_codes)
lengths = Array(Int64,n_codes+1)
for i = 1:n_codes
  fname = string(parsed_args["prefix"], "_", code_map[i,1],
                 parsed_args["suffix"])
  code_files[i] = open(fname,"w")
  code_counts[i] = 0
  lengths[i] = length(code_map[i,2])
end
if KEEP_UNMATCHED
  fname = string(parsed_args["prefix"], "_unmatched", parsed_args["suffix"])
  code_files[n_codes+1] = open(fname,"w")
end
lengths[n_codes+1] = 0

function match_barcode(sequence,barcodes,lengths,mismatches,partial,beginning)
  n_barcodes = length(barcodes)
  misses = Array(Int64,n_barcodes,partial+1)

  if beginning
    for p = 0 : partial
      for code = 1 : n_barcodes
        misses[code,p+1] = p
        k = lengths[code]
        for i = 1 : k - p
          misses[code,p+1] += barcodes[code][i+p] != sequence[i]
        end
        if misses[code,p+1] == 0
          return code, p, misses[code,p+1]
        end
      end
    end
  else
    for p = 0 : partial
      for code = 1 : n_barcodes
        misses[code,p+1] = p
        k = lengths[code]
        for i = 1 : k - p
          misses[code,p+1] += barcodes[code][i] != sequence[end-k+i+p]
        end
        if misses[code,p+1] == 0
          return code, p, misses[code,p+1]
        end
      end
    end
  end

  mincode,minpart = findn(misses .== minimum(misses))
  if misses[mincode[1],minpart[1]] > mismatches
    result = 0
  else
    result = mincode[1]
  end
  return result, minpart[1]-1, misses[mincode[1],minpart[1]]
end

overhang_barcodes = ["TAACAG", "TAACA", "TAAC", "TAA"]
overhang_lengths = [6, 5, 4, 3]

overhang_fails = 0
unmatched_count = 0

tic()
while !eof(input)
  name1 = chomp(readline(input))
  seq = chomp(readline(input))
  name2 = chomp(readline(input))
  quality = chomp(readline(input))

  if length(seq) < 51
    continue
  end

  # trim off first 9 (HiSeq initiation) and last 17 (adapter) bases
  seq = seq[10:end-17]
  quality = quality[10:end-17]

  trimcode,part,miss = match_barcode(seq,overhang_barcodes,overhang_lengths,
                                     1,1,false)
  if trimcode > 0
    seq = seq[1:20-overhang_lengths[trimcode]+7]
    quality = quality[1:20-overhang_lengths[trimcode]+7]
  else
    # don't output this read
    overhang_fails += 1
    continue
  end

  idx,part,misc = match_barcode(seq,code_map[:,2],lengths,1,1,true)

  if idx == 0
    unmatched_count += 1
    if KEEP_UNMATCHED
      # assign to "unmatched" file
      idx = n_codes + 1
    else
      # skip this read
      continue
    end
  else
    code_counts[idx] += 1
  end
  @printf(code_files[idx],"%s\n",name1)
  @printf(code_files[idx],"%s\n",seq[lengths[idx]+1:end])
  @printf(code_files[idx],"%s\n",name2)
  @printf(code_files[idx],"%s\n",quality[lengths[idx]+1:end])
end
elapsed_seconds = toq()

for i = 1:n_codes
  close(code_files[i])
end
if KEEP_UNMATCHED
  close(code_files[n_codes+1])
end

total_reads = sum(code_counts) + unmatched_count + overhang_fails
decoded_count = sum(code_counts)
@printf("pool\treads\tpercent_of_total\tpercent_of_decoded\n")
for i = 1 : n_codes
  @printf("%s\t%i\t%.3f\t%.3f\n",
          code_map[i,1],code_counts[i],code_counts[i]/total_reads*100,
          code_counts[i]/decoded_count*100)
end
@printf("%s\t%i\t%.3f\t%.3f\n",
        "unmatched",unmatched_count,unmatched_count/total_reads*100,
        unmatched_count/decoded_count*100)
@printf("%s\t%i\t%.3f\t%.3f\n",
        "removed",overhang_fails,overhang_fails/total_reads*100,
        overhang_fails/decoded_count*100)

@printf("#time %f seconds\n",elapsed_seconds)
