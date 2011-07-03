INDEX_MAP = {
    "ATCACG"  =>  "AR01",
    "CGATGT"  =>  "AR02",
    "TTAGGC"  =>  "AR03",
    "TGACCA"  =>  "AR04",
    "ACAGTG"  =>  "AR05",
    "GCCAAT"  =>  "AR06",
    "CAGATC"  =>  "AR07",
    "ACTTGA"  =>  "AR08",
    "GATCAG"  =>  "AR09",
    "TAGCTT"  =>  "AR10",
    "GGCTAC"  =>  "AR11",
    "CTTGTA"  =>  "AR12"
  }

read_1_file = ARGV[0]
index_file  = ARGV[1]
indices     = ARGV[2,ARGV.length-1]

basename    = read_1_file.gsub("_qseq.txt", "")

INDEX_FILES = {}

File.open(read_1_file, "r") do |read_1_in|
  File.open(index_file, "r") do |index_in|
    while read_1_line = read_1_in.gets
      index_line = index_in.gets.chomp
      index = index_line.split("\t")[8][0,6]
      next unless INDEX_MAP.has_key? index
      index_name = INDEX_MAP[index]
      next unless indices.include? index_name
      unless INDEX_FILES.has_key? index
        INDEX_FILES[index] = File.open(basename + "_" + index_name + "_qseq.txt","w")
      end
      INDEX_FILES[index].puts read_1_line
    end
  end
end

for key in INDEX_FILES.keys
  INDEX_FILES[key].close
end