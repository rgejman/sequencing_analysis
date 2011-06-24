INDEX_MAP = {
    "CGATGT"  =>  "AR02",
    "TGACCA"  =>  "AR04",
    "ACAGTG"  =>  "AR05",
    "GCCAAT"  =>  "AR06",
    "CAGATC"  =>  "AR07"
    "CTTGTA"  =>  "AR12"
  }

read_1_file = ARGV[0]
index_file  = ARGV[1]

basename    = read_1_file.gsub(".txt", "")

INDEX_FILES = {}

File.open(read_1_file, "r") do |read_1_in|
  File.open(index_file, "r") do |index_in|
    while read_1_line = read_1_in.gets
      index_line = index_in.gets.chomp
      index = index_line.split("\t")[8]
      next unless INDEX_MAP.has_key? index
      unless INDEX_FILES.has_key? index
        INDEX_FILES[index] = File.open(basename + "_" + INDEX_MAP[index] + ".txt","w")
      end
      INDEX_FILES[index].puts read_1_line
    end
  end
end

for key in INDEX_FILES.keys
  INDEX_FILES[key].close
end