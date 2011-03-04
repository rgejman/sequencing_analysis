#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 12

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM rna_seq_alignment ORDER BY created_at desc")
res.each_hash do |rna_seq_alignment|
  person                        = rna_seq_alignment["person"]
  output_folder_name            = person + "_" + rna_seq_alignment["sample"]
  running_file                  = running_file(output_folder_name, "cufflinks")
  output_folder_path            = "#{TOPHAT_FOLDER}/#{person}/#{output_folder_name}"
  accepted_hits                 = "#{output_folder_path}/accepted_hits.bam"
  junctions                     = "#{output_folder_path}/junctions.bed"
  transcripts                   = "#{output_folder_path}/transcripts"
  puts "Checking for existence of #{output_folder_path}"
  puts "Checking for existence of #{accepted_hits}"
  puts "Checking for existence of #{junctions}"
  next if File.exists? "#{output_folder_path}/transcripts/transcripts.gtf"
  next if File.exists? running_file
  next unless File.exists? output_folder_path #Tophat has not yet run.
  next unless File.exists? accepted_hits
  next unless File.exists? junctions

  REF_TRANSCRIPTS_FILE = "#{USEFUL_BED_FILES}/Mus_musculus.NCBIM37.61.for-tophat.gtf"
  `touch #{running_file}`
  begin
    cmd = "cufflinks -G #{REF_TRANSCRIPTS_FILE} -p #{NUM_THREADS} -o #{transcripts} #{accepted_hits}" # -r #{BOWTIE_INDEXES}/#{GENOME}.fa
    puts cmd
    `#{cmd}`
    Dir.chdir("#{output_folder_path}/transcripts")
    `cuffcompare -r #{REF_TRANSCRIPTS_FILE} -o cuffcompare transcripts.gtf`
    
    # Replace cufflinks gene_ids with the gene symbols
    gene_expr = File.readlines("#{transcripts}/genes.expr").collect{|l| l.chomp.split("\t")}
    loci = {}
    File.readlines("#{transcripts}/cuffcompare.loci").each{|l| t=l.split("\t"); loci[t[0]]=t[2].split("|").first}
    File.open("#{transcripts}/genes.symbols.expr","w") do |f|
      f.puts gene_expr.shift.join("\t")
      for line in gene_expr
        line[0] = loci[line[0]]
        f.puts line.join("\t")
      end
    end
    
  rescue => e
    throw e
  ensure
    FileUtils.rm("#{output_folder_path}/transcripts",    :force=>true)
    FileUtils.rm(running_file,                    :force=>true)
    conn.close
  end
  break
end