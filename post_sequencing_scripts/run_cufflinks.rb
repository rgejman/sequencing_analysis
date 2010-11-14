#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 12

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM rna_seq_alignment ORDER BY created_at desc")
res.each_hash do |rna_seq_alignment|
  output_folder_name            = rna_seq_alignment["person"] + "_" + rna_seq_alignment["sample"]
  running_file                  = running_file(output_folder_name, "cufflinks")
  output_folder_path            = "#{TOPHAT_FOLDER}/#{output_folder_name}"
  accepted_hits                 = "#{output_folder_path}/accepted_hits.bam"
  junctions                     = "#{output_folder_path}/junctions.bed"
  puts "Checking for existence of #{output_folder_path}"
  puts "Checking for existence of #{accepted_hits}"
  puts "Checking for existence of #{junctions}"
  next if File.exists? "#{output_folder_path}/transcripts/transcripts.gtf"
  next if File.exists? running_file
  next unless File.exists? output_folder_path #Tophat has not yet run.
  next unless File.exists? accepted_hits
  next unless File.exists? junctions

  REF_TRANSCRIPTS_FILE = "#{USEFUL_BED_FILES}/mm9.ucsc.genes.gtf"
  `touch #{running_file}`
  begin
    cmd = "cufflinks -p #{NUM_THREADS} -r #{BOWTIE_INDEXES}/#{GENOME}.fa -o #{output_folder_path}/transcripts #{accepted_hits}"
    puts cmd
    `#{cmd}`
    `cd transcripts`
    `cuffcompare -r #{REF_TRANSCRIPTS_FILE} -o cuffcompare transcripts.gtf`
  rescue => e
    throw e
  ensure
    FileUtils.rm("#{output_folder_path}/transcripts",    :force=>true)
    FileUtils.rm(running_file,                    :force=>true)
    conn.close
  end
  break
end