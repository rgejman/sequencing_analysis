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
  
for key in INDEX_MAP.keys
  a = key.split("")
  for i in 0...a.length
    index = a[0,i].to_s + "." + a[i+1,6].to_s
    INDEX_MAP[index] = INDEX_MAP[index]
    puts index
  end
end

paired      = ARGV[0].to_i == 2
read_1_file = ARGV[1]

if !paired
  index_file  = ARGV[2]
  indices     = ARGV[3,ARGV.length-1]
  basename    = read_1_file.gsub("_qseq.txt", "")
  
else
  read_2_file = ARGV[2]
  index_file  = ARGV[3]
  indices     = ARGV[4,ARGV.length-1]
  basename    = read_1_file.gsub("_1_qseq.txt", "")
end


INDEX_FILES = {}

read_1_in = File.open(read_1_file, "r")
read_2_in = File.open(read_2_file, "r") if paired

File.open(index_file, "r") do |index_in|
  while read_1_line = read_1_in.gets
    index_line = index_in.gets.chomp
    read_2_line = read_2_in.gets if paired
    index = index_line.split("\t")[8][0,6]
    next unless INDEX_MAP.has_key? index
    index_name = INDEX_MAP[index]
    next unless indices.include? index_name
    unless INDEX_FILES.has_key? index + "_read_1"
      if !paired
        INDEX_FILES[index + "_read_1"] = File.open(basename + "_" + index_name + "_qseq.txt","w")
      else
        INDEX_FILES[index + "_read_1"] = File.open(basename + "_" + index_name + "_1_qseq.txt","w")
        INDEX_FILES[index + "_read_2"] = File.open(basename + "_" + index_name + "_2_qseq.txt","w")
      end
    end
    INDEX_FILES[index + "_read_1"].puts read_1_line
    INDEX_FILES[index + "_read_2"].puts read_2_line if paired
    
  end
end

read_1_in.close
read_2_in.close if paired

for key in INDEX_FILES.keys
  INDEX_FILES[key].close
end