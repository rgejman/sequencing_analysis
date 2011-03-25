#!/usr/bin/env ruby -wKU
## For each gene: calculates "Pausing Index" as defined by Muse et al (2007 Nat Genetics 29, 1507-1511)
## Pausing_index = log2((PolII at TSS +/- 250bp) - (PolII +500bp to end of gene))
require 'wigreader'
require 'pp'

class << Math
  def log2(n); log(n) / log(2); end
end

WIG_FILE            = ARGV[0] ## This should be a RNA Pol II chip-seq wiggle file
US_TSS_OFFSET       = 250 # +/- this number is the area defined to be the paused-pol II area.
TSS_DS_OFFSET       = 500 # TSS + this number until the end of the gene is the "DS" area

GENE_BED_FILE       = "/home/tarakhovsky/genomics/useful_bed_files/mm9.refseq.genes.uniq_longer.bed"
OUTPUT_FILE         = File.basename(WIG_FILE).split(".")[0] + ".pausing_index.txt" #e.g. "Eugene_WT_H3K4me3_CD4.sorted.bam.wig" to "Eugene_WT_H3K4me3_CD4.profile.6000.80.txt"

block_len = block_len.to_i
genes = File.readlines(GENE_BED_FILE).collect{|l| g={}; g[:chr],g[:start],g[:end],g[:symbol],tmp,g[:strand]= l.chomp.split("\t");g}
genes = genes.reject {|g| g[:chr] =~ /random/ }
reader = WigReader.new(WIG_FILE)

for gene in genes
  gene[:end]    = gene[:end].to_i
  gene[:start]  = gene[:start].to_i
  pausing_start = gene[:start] - US_TSS_OFFSET
  pausing_end   = gene[:start] + TSS_DS_OFFSET
  ds_start      = gene[:start] + TSS_DS_OFFSET + 1
  ds_end        = gene[:end]
  if gene[:strand] == "-"
    pausing_start = gene[:start] + US_TSS_OFFSET
    pausing_end   = gene[:start] - TSS_DS_OFFSET
    ds_start      = gene[:start] - TSS_DS_OFFSET - 1
    ds_end        = gene[:end]
  end
  begin
    gene[:pausing_index]      = 0 # this means it could not be calculated.
    gene[:pausing_area_rpkm]  = reader.fpkm(gene[:chr], pausing_start, pausing_end)
    gene[:ds_area_rpkm]       = reader.fpkm(gene[:chr], ds_start, ds_end)
    diff = gene[:pausing_area_rpkm] - gene[:ds_area_rpkm]
    if diff <= 0
      gene[:pausing_index]    = 0
    else
      gene[:pausing_index]    = Math.log2(diff)
    end
  rescue => e
    puts e
    puts "#{gene[:symbol]} could not be found in the wig file."
  end
end
File.open(OUTPUT_FILE, "w") do |f|
  for gene in genes
    f.puts "#{gene[:symbol]}\t#{gene[:chr]}\t#{gene[:start]}\t#{gene[:end]}\t#{gene[:strand]}\t#{gene[:pausing_index]}\t#{gene[:pausing_area_rpkm]}\t#{gene[:ds_area_rpkm]}"
  end
end