#!/usr/bin/env ruby -wKU
## For each gene: calculates "Pausing Index" as defined by Muse et al (2007 Nat Genetics 29, 1507-1511)
## Pausing_index = (PolII at TSS +/- 250bp) / (PolII +500bp to end of gene)
require 'wigreader'
require 'pp'

def remove_smaller_genes(arr)
  arr.uniq! # To remove genes that have the same TSS and TES
  f_arr = arr.clone
  for a in (0...arr.length)
    for b in ((a+1)...arr.length)
      next unless a[:gene] == b[:gene]
      a_len = (g[:start] - g[:end]).abs
      b_len = (g[:start] - g[:end]).abs
      if a_len < b_len
        f_arr.remove(a)
      elsif
        f_arr.remove(b)
      end
    end
  end
  return f_arr
end

WIG_FILE            = ARGV[0] ## This should be a RNA Pol II chip-seq wiggle file
AROUND_TSS_AREA     = 250 # +/- this number is the area defined to be the paused-pol II area.
TSS_DS_OFFSET       = 500 # TSS + this number until the end of the gene is the "DS" area

GENE_BED_FILE       = "/home/tarakhovsky/genomics/useful_bed_files/mm9.refseq.genes.gtf"
output_file         = File.basename(WIG_FILE).split(".")[0] + ".pausing_index.txt" #e.g. "Eugene_WT_H3K4me3_CD4.sorted.bam.wig" to "Eugene_WT_H3K4me3_CD4.profile.6000.80.txt"

block_len = block_len.to_i
genes = File.readlines(GENE_BED_FILE).collect{|l| g={}; g[:chr],g[:start],g[:end],g[:gene],tmp,g[:strand]= l.chomp.split("\t");g}
genes = genes.reject {|g| g[:chr] =~ /random/ }

## Some genes have multiple entries in this file because of differing splice patterns.
## Choose the longer genes (as defined by TSS and Transcription termination site)
## Why? Because you're less likely to have false positive paused polymerase genes 
genes = remove_smaller_genes(genes)

reader = WigReader.new(WIG_FILE)

for gene in genes
  gene[:end]    = gene[:end].to_i
  gene[:start]  = gene[:start].to_i
  pausing_start = gene[:start] - AROUND_TSS_AREA
  pausing_end   = gene[:start] + AROUND_TSS_AREA
  ds_start      = gene[:start] + TSS_DS_OFFSET
  ds_end        = gene[:end]
  if gene[:strand] == "-"
    pausing_start = gene[:start] + AROUND_TSS_AREA
    pausing_end   = gene[:start] - AROUND_TSS_AREA
    ds_start      = gene[:start] - TSS_DS_OFFSET
    ds_end        = gene[:end]
  end
  begin
    gene[:pausing_area_rpkm]  = reader.fpkm(gene[:chr], pausing_start, pausing_end)
    gene[:ds_area_rpkm]       = reader.fpkm(gene[:chr], ds_start, ds_end)
    gene[:pausing_index]      = gene[:pausing_area_rpkm] / gene[:ds_area_rpkm]
  rescue => e
    puts e
    puts "#{gene[:gene]} could not be found in the wig file."
  end
end

File.open(output_file, "w") do |f|
  f.print "symbol\tchr\tstart\tend\tstrand\tpausing_index\tpromotor_rpkm\tds_rpkm"
  f.print "\n"
  for gene in genes
    f.puts "#{gene[:gene]}\t#{gene[:chr]}\t#{gene[:start]}\t#{gene[:end]}\t#{gene[:strand]}\t#{gene[:pausing_index]}\t#{gene[:pausing_area_rpkm]}\t#{gene[:ds_area_rpkm]}"
  end
end