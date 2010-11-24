#!/usr/bin/env ruby -wKU
require 'wigreader'
require 'pp'

WIG_FILE            = ARGV[0]
GENE_BED_FILE       = ARGV[1]
DISTANCE_AROUND_TSS = ARGV[2].to_i
BLOCKS              = ARGV[3].to_i
OUTPUT_FILE         = ARGV[4]

block_len = DISTANCE_AROUND_TSS.to_f/BLOCKS
raise "'Distance around TSS' must be evenly divisible by 'blocks'" if block_len != block_len.round
block_len = block_len.to_i
genes = File.readlines(GENE_BED_FILE).collect{|l| g={}; g[:chr],g[:start],g[:end],g[:gene],tmp,g[:strand]= l.chomp.split("\t");g}
genes = genes.reject {|g| g[:chr] =~ /random/ }
reader = WigReader.new(WIG_FILE)

for gene in genes
  gene[:end]    = gene[:end].to_i
  gene[:start]  = gene[:start].to_i
  d             = gene[:end] - gene[:start]
  to_add        = (DISTANCE_AROUND_TSS - d)/2
  gene[:start]  -= to_add
  gene[:end]    += to_add
  gene[:blocks] = []
  begin
    for i in (0...BLOCKS)
      offset = i*block_len
      if gene[:strand] == "+"
        start   = gene[:start]+offset
        last    = start+block_len-1
      else
        last    = gene[:end]-offset
        start   = last-block_len+1
      end
      #puts "Will read FPKM for: #{gene[:gene]} #{gene[:strand]} #{gene[:chr]} #{start}-#{last} #{last-start}"
      gene[:blocks][i] = reader.fpkm(gene[:chr], start, last)
    end
  rescue => e
    puts e
    puts "#{gene[:gene]} could not be found in the wig file. Try adjusting the distance_around_tss"
    for i in (0...BLOCKS)
      gene[:blocks][i] = 0.0
    end
  end
end

File.open(OUTPUT_FILE, "w") do |f|
  f.print "symbol\tchr\tstart\tend\tstrand"
  for i in 0...genes[0][:blocks].length
    f.print "\tb#{i}"
  end
  for gene in genes
    f.print "#{gene[:gene]}\t#{gene[:chr]}\t#{gene[:start]}\t#{gene[:end]}\t#{gene[:strand]}"
    for block in gene[:blocks]
      f.print "\t#{block}"
    end
    f.print "\n"
  end
end