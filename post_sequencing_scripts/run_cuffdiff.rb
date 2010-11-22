#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 24

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM rna_seq_analysis_pairs ORDER BY created_at desc")
res.each_hash do |pair|
  foreground_id     = pair["foreground_rna_seq_alignment_id"]
  background_id     = pair["background_rna_seq_alignment_id"]
  foreground        = conn.query("SELECT * FROM rna_seq_alignment WHERE id = '#{foreground_id}'").fetch_hash
  background        = conn.query("SELECT * FROM rna_seq_alignment WHERE id = '#{background_id}'").fetch_hash
  
  foreground_name   = foreground["person"] + "_" + foreground["sample"]
  background_name   = background["person"] + "_" + background["sample"]
  
  foreground_folder = "#{TOPHAT_FOLDER}/#{foreground_name}"
  background_folder = "#{TOPHAT_FOLDER}/#{background_name}"
  
  output_folder_name  = "#{foreground['person']}_#{foreground['sample']}_#{background['sample']}"
  output_folder_path  = "#{DIFF_EXPR_FOLDER}/#{output_folder_name}"
  
  labels              = background["sample"] + "," + foreground["sample"]
  
  next if File.exists? output_folder_path
  
  print "Checking existence of #{foreground_folder}/transcripts/transcripts.gtf ... "
  next unless File.exists? "#{foreground_folder}/transcripts/transcripts.gtf"
  puts "exists"
  print "Checking existence of #{background_folder}/transcripts/transcripts.gtf ... "
  next unless File.exists? "#{background_folder}/transcripts/transcripts.gtf"
  print "exists"
  running_file      = running_file(output_folder_name, "cuffdiff")
  
  REF_TRANSCRIPTS_FILE = "#{USEFUL_BED_FILES}/mm9.ensemble.genes.for.cuffdiff.fixed.gtf"
  `touch #{running_file}`
  begin
    cmd = "cuffdiff -o #{output_folder_path} -p #{NUM_THREADS} -L #{labels} -r #{BOWTIE_INDEXES}/#{GENOME}.fa #{REF_TRANSCRIPTS_FILE} #{background_folder}/accepted_hits.bam #{foreground_folder}/accepted_hits.bam"
    puts cmd
    `#{cmd}`
  rescue => e
    throw e
  ensure
    FileUtils.rm(running_file,                    :force=>true)
    FileUtils.rm(output_folder_path,              :force=>true)
    conn.close
  end
  break
end