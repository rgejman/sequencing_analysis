#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 12

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM rna_seq_alignment ORDER BY created_at desc")
res.each_hash do |rna_seq_alignment|
  person =  rna_seq_alignment["person"]
  output_folder_name  = person + "_" + rna_seq_alignment["sample"]
  output_folder_path  = "#{TOPHAT_FOLDER}/#{person}/#{output_folder_name}"
  next if File.exists? output_folder_path #this has already been analyzed. 
  files_res = conn.query("SELECT * FROM rna_seq_pairs WHERE rna_seq_alignment_id = #{rna_seq_alignment['id']}")
  reads = []
  # for each file that comprises this RNA-Seq run (there may be multiple)
  files_res.each_hash do |file|
    if file["paired"].to_i == 1
      f1 = "#{FASTQ_RNA_SEQ_FOLDER}/#{person}/#{file['name']}.1.txt"
      f2 = "#{FASTQ_RNA_SEQ_FOLDER}/#{person}/#{file['name']}.2.txt"
      reads << [f1,f2]
    else
      reads << "#{FASTQ_RNA_SEQ_FOLDER}/#{person}/#{file['name']}.txt"
    end
  end
  read_length           = rna_seq_alignment["read_length"].to_i
  mean_fragment_length  = rna_seq_alignment["mean_fragment_length"].to_i
  genome                = rna_seq_alignment["genome"]
  genome                = GENOMES[genome]
  
  mean_dist_arg         = ""
  if reads[0].length == 2 #the first "reads entry" has 2 files, so it's paired.
    files_arg = reads.collect {|p| p[0] }.join(",") + " " + reads.collect {|p| p[1] }.join(",")
    mean_dist           = mean_fragment_length - (read_length*2)# - 70 #70 accounts for the illumina primers
    mean_dist_arg       = "-r #{mean_dist}" #-r = mean distance between ends of paired reads.
  else
    files_arg = reads.collect {|p| p[0] }.join(",")
  end
  GTF_FILE_ARG = "-G #{USEFUL_BED_FILES}/rna_seq/Mus_musculus.NCBIM37.61.for-tophat.gtf"
  LIBRARY_TYPE_ARG = "--library-type fr-unstranded"
  OUTPUT_FOLDER_ARG = "-o #{output_folder_path}"
  begin
    cmd = "tophat -a 7 -p #{NUM_THREADS} #{mean_dist_arg} #{GTF_FILE_ARG} #{LIBRARY_TYPE_ARG} #{OUTPUT_FOLDER_ARG} #{genome} #{files_arg}"
    puts cmd
    `#{cmd}`
    `bamtools index -in #{output_folder_path}/accepted_hits.bam`
  rescue => e
    throw e
  ensure
    FileUtils.rm(output_folder_path,    :force=>true)
    conn.close
  end
  break
end