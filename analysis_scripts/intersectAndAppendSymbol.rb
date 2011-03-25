#!/usr/bin/env ruby -wKU

## Take in 2 files
## 1: A file with loci information
## 2: a BED file with gene information
## Output a file with the loci_file information + the symbol from the BED file.
## Except that if multiple lines in the Loci file match a single entry in the BED file, only output one of the peaks.
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
ARGS_REQUIRED = 5
unless ARGV.length == ARGS_REQUIRED
  puts "Usage:\t\tintersectAndAppendSymbol.rb loci_file chr_col start_col end_col bed_file"
  files = Dir.glob("#{USEFUL_BED_FILES}/*.bed").collect {|f| f.split("/").last}
  puts ""
  puts "Region files:\t#{files.join("\n\t\t")}"
  exit
end

used_symbols = []

LOCI_FILE       = ARGV[0]
CHR_COL         = ARGV[1].to_i #the column which has the chromosome in the loci file
START_COL       = ARGV[2].to_i #the column which has the start position for the loci file
END_COL         = ARGV[3].to_i #the column which has the end position for the loci file
BED_FILE        = ARGV[4]

BED             = {}
bed_lines       = File.readlines(BED_FILE).collect{|l| t=l.chomp.split("\t"); BED[t[0]] ||= []; BED[t[0]] << t }

loci_lines      = File.readlines(LOCI_FILE).collect{|l| l.chomp.split("\t")}
puts loci_lines.shift.join("\t") + "\tSymbol" #header
for loci in loci_lines
  chr       = loci[CHR_COL]
  start_pos = loci[START_COL].to_i
  end_pos   = loci[END_COL].to_i
  
  ## POSSIBLE BUG HERE: BED files usually have coordinates from lower --> higher, regardless of strand.
  
  if start_pos > end_pos
    tmp       = end_pos
    end_pos   = start_pos
    start_pos = tmp
  end
  for gene in BED[chr]
    s = gene[1].to_i
    e = gene[2].to_i
    if (s >= start_pos and s <= end_pos) or (e >= start_pos and e <= end_pos) or (s <= start_pos and e >= end_pos)
      # intersect!
      next if used_symbols.include? gene[3]
      puts loci.join("\t") + "\t" + gene[3] #append the gene symbol
      used_symbols << gene[3]
      break
    end
  end
end