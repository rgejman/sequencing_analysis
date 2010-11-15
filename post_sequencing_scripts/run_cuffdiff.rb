#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 12

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
  
  labels            = foreground["sample"] + "," + background["sample"]
  
  print "Checking existence of #{foreground_folder}/transcripts.gtf ... "
  next unless File.exists? "#{foreground_folder}/transcripts.gtf"
  puts "exists"
  print "Checking existence of #{background_folder}/transcripts.gtf ... "
  next unless File.exists? "#{background_folder}/transcripts.gtf"
  print "exists"
  running_file      = running_file(output_folder_name, "cuffdiff")
  
  REF_TRANSCRIPTS_FILE = "#{USEFUL_BED_FILES}/mm9.ucsc.genes.gtf"
  `touch #{running_file}`
  begin
    cmd = "cuffdiff -p #{NUM_THREADS} -L #{labels} -r #{BOWTIE_INDEXES}/#{GENOME}.fa #{REF_TRANSCRIPTS_FILE} #{foreground_folder}/accepted_hits.bam #{background_folder}/accepted_hits.bam"
    puts cmd
    `#{cmd}`
  rescue => e
    throw e
  ensure
    FileUtils.rm("#{output_folder_path}/transcripts",    :force=>true)
    FileUtils.rm(running_file,                    :force=>true)
    conn.close
  end
  break
end