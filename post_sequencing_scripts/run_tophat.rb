#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 12

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM rna_seq_alignment ORDER BY created_at desc")
res.each_hash do |rna_seq_alignment|
  output_folder_name  = rna_seq_alignment["person"] + "_" + rna_seq_alignment["sample"]
  output_folder_path  = "#{TOPHAT_FOLDER}/#{output_folder_name}"
  next if File.exists? output_folder_path #this has already been analyzed. 
  files_res = conn.query("SELECT * FROM rna_seq_pairs WHERE rna_seq_alignment_id = #{rna_seq_alignment['id']}")
  reads = []
  # for each file that comprises this RNA-Seq run (there may be multiple)
  files_res.each_hash do |file|
    if file["paired"].to_i == 1
      f1 = "#{FASTQ_RNA_SEQ_FOLDER}/#{file['name']}.1.txt"
      f2 = "#{FASTQ_RNA_SEQ_FOLDER}/#{file['name']}.2.txt"
      reads << [f1,f2]
    else
      reads << "#{FASTQ_RNA_SEQ_FOLDER}/#{file['name']}.txt"
    end
  end
  read_length           = rna_seq_alignment["read_length"].to_i
  mean_fragment_length  = rna_seq_alignment["mean_fragment_length"].to_i
  mean_dist_arg         = ""
  if reads[0].length == 2 #the first "reads entry" has 2 files, so it's paired.
    files_arg = reads.collect {|p| p[0] }.join(",") + " " + reads.collect {|p| p[1] }.join(",")
    mean_dist           = mean_fragment_length - (read_length*2) - 70 #70 accounts for the illumina primers
    mean_dist_arg       = "-r #{mean_dist}" #-r = mean distance between ends of paired reads.
  else
    files_arg = reads.collect {|p| p[0] }.join(",")
  end
  GTF_FILE_ARG = "-G #{USEFUL_BED_FILES}/mm9.ucsc.genes.gtf"
  LIBRARY_TYPE_ARG = "--library-type fr-unstranded"
  OUTPUT_FOLDER_ARG = "-o #{output_folder_path}"
  begin
    cmd = "tophat -p #{NUM_THREADS} #{mean_dist_arg} #{GTF_FILE_ARG} #{LIBRARY_TYPE_ARG} #{OUTPUT_FOLDER_ARG} #{GENOME} #{files_arg}"
    puts cmd
    `#{cmd}`
  rescue => e
    throw e
  ensure
    FileUtils.rm(output_folder_path,    :force=>true)
    conn.close
  end
  break
end