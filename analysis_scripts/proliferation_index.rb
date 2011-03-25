## Calculate the Proliferation Indices of a sample
## There are two: Positive Proliferation Index (PPI)
##                Negative Proliferation Index (NPI)

# XPI = sum(proliferation_genes.fpkm)/proliferation_genes.length

# Proliferation Index (PI) = PPI/NPI = The ratio of positive proliferation genes to negative proliferation genes
# The higher the PI, the more the cells are proliferating (can we show this?)

## Which are pro and which are anti proliferation genes?
## http://www.geneontology.org has two GO categories for this:
## In mus musculus:
## 1. Positive Regulation of Cell Proliferation (GO:0008285) => 495 genes/proteins
## 2. Negative Regulation of Cell Proliferation (GO:0008284) => 345 genes/proteins

# INPUTS:
# 1: SAMPLE_FILE
# The file that lists FPKM values per gene.
# This should be a tab delimited file with a gene-symbol or ID column and a FPKM column (with headers)
# e.g.:
# Symbol  FPKM
# Zeb1  15.347
# Dcp2  16.5366
# 2: A list of positive proliferation gene symbols/IDS
# 3: A list of negative proliferation gene symbols/IDS

# The gene symbols/IDs of the PPG and NPG files must match the gene symbol/IDs of the FPKM file

SAMPLE_FILE     = ARGV[0]
GENE_ID_COLUMN  = ARGV[1].to_i #the column which has the gene ID/symbol
FPKM_COLUMN     = ARGV[2].to_i #the column which has the FPKM value
PPG             = File.readlines(ARGV[3]).collect{|l| l.chomp} # Positive Proliferation Genes
NPG             = File.readlines(ARGV[4]).collect{|l| l.chomp} # Negative Proliferation Genes

fpkms         = File.readlines(SAMPLE_FILE).collect {|l| s=l.split("\t"); {:id=>s[GENE_ID_COLUMN],:fpkm=>s[FPKM_COLUMN].to_f}}

pp_fpkm = 0
np_fpkm = 0

for gene in fpkms
  if PPG.include? gene[:id]
    pp_fpkm += gene[:fpkm]
  end  
  if NPG.include? gene[:id]
    np_fpkm += gene[:fpkm]
  end
end

ppi = pp_fpkm / PPG.length
npi = np_fpkm / NPG.length

puts "PP_FPKM\t#{pp_fpkm}"
puts "NP_FPKM\t#{np_fpkm}"
puts ""

puts "PPG\t#{PPG.length}"
puts "NPG\t#{NPG.length}"
puts ""

puts "PPI\t#{ppi}"
puts "NPI\t#{npi}"
puts ""

puts "PI\t#{ppi/npi}"