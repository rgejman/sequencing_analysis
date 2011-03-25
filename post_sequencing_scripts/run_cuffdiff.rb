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
  
  f_person          = foreground["person"]
  b_person          = background["person"]
  
  foreground_name   = f_person + "_" + foreground["sample"]
  background_name   = b_person + "_" + background["sample"]
  
  foreground_folder = "#{TOPHAT_FOLDER}/#{f_person}/#{foreground_name}"
  background_folder = "#{TOPHAT_FOLDER}/#{b_person}/#{background_name}"
  
  output_folder_name  = "#{f_person}_#{foreground['sample']}_#{background['sample']}"
  output_folder_path  = "#{DIFF_EXPR_FOLDER}/#{f_person}/#{output_folder_name}"
  
  labels              = background["sample"] + "," + foreground["sample"]
  
  next if File.exists? output_folder_path
  
  print "Checking existence of #{foreground_folder}/transcripts/transcripts.gtf ... "
  next unless File.exists? "#{foreground_folder}/transcripts/transcripts.gtf"
  puts "exists"
  print "Checking existence of #{background_folder}/transcripts/transcripts.gtf ... "
  next unless File.exists? "#{background_folder}/transcripts/transcripts.gtf"
  puts "exists"
  running_file      = running_file(output_folder_name, "cuffdiff")
  
  #-r #{BOWTIE_INDEXES}/#{GENOME}.fa
  REF_TRANSCRIPTS_FILE = "#{USEFUL_BED_FILES}/rna_seq/Mus_musculus.NCBIM37.61.for-tophat.gtf"
  `touch #{running_file}`
  begin
    `mkdir -p #{output_folder_path}`
    cmd = "cuffdiff -o #{output_folder_path} -p #{NUM_THREADS} -L #{labels} #{REF_TRANSCRIPTS_FILE} #{background_folder}/accepted_hits.bam #{foreground_folder}/accepted_hits.bam"
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