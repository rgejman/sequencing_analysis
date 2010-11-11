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
  tophat_output_folder_path     = "#{TOPHAT_FOLDER}/#{output_folder_name}"
  cufflinks_output_folder_path  = "#{CUFFLINKS_FOLDER}/#{output_folder_name}"
  accepted_hits                 = "#{tophat_output_folder_path}/accepted_hits.bam"
  junctions                     = "#{tophat_output_folder_path}/junctions.bed"
  puts "Checking for existence of #{tophat_output_folder_path}"
  puts "Checking for existence of #{accepted_hits}"
  puts "Checking for existence of #{junctions}"
  next unless File.exists? tophat_output_folder_path #Tophat has not yet run.
  next unless File.exists? accepted_hits
  next unless File.exists? junctions
  files_res = conn.query("SELECT * FROM rna_seq_pairs WHERE rna_seq_alignment_id = #{rna_seq_alignment['id']}")
  paired = files_res.fetch_hash["paired"] == "1"
  read_length           = rna_seq_alignment["read_length"].to_i
  mean_fragment_length  = rna_seq_alignment["mean_fragment_length"].to_i
  mean_dist_arg         = ""
  if paired #the first "reads entry" has 2 files, so it's paired.
    mean_dist           = mean_fragment_length - (read_length*2) - 70 #70 accounts for the illumina primers
    mean_dist_arg       = "-m #{mean_dist}" #-r = mean distance between ends of paired reads.
  end
  
  #GTF_FILE_ARG = "-G #{USEFUL_BED_FILES}/mm9.ucsc.genes.gtf"
  #LIBRARY_TYPE_ARG = "--library-type fr-unstranded" # this is the default
  `touch #{running_file}`
  begin
    cmd = "tophat -p #{NUM_THREADS} #{mean_dist_arg} -o #{cufflinks_output_folder_path} -r #{BOWTIE_INDEXES}/#{GENOME}.fa #{accepted_hits}"
    puts cmd
    `#{cmd}`
  rescue => e
    throw e
  ensure
    FileUtils.rm(cufflinks_output_folder_path,    :force=>true)
    FileUtils.rm(running_file,                    :force=>true)
    conn.close
  end
  break
end