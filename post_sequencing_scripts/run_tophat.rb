#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 16

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM rna_seq_alignment ORDER BY created_at desc")
res.each_hash do |rna_seq_alignment|
  output_folder_name  = rna_seq_alignment["person"] + "_" + rna_seq_alignment["sample"]
  output_folder_path  = "#{TOPHAT_FOLDER}/#{analysis_folder_name}"
  next if File.exists? output_folder_path #this has already been analyzed. 
  files_res = conn.query("SELECT * FROM rna_seq_pairs WHERE rna_seq_alignment_id = #{rna_seq_alignment['id']}")
  pairs = []
  # for each file that comprises this RNA-Seq run (there may be multiple)
  files_res.each_hash do |file|
    if file["paired"].to_i == 1
      pairs << [file["name"] + ".1.txt", file["name"] + ".2.txt"]
    else
      pairs << file["name"] + ".txt"
    end
  end
  read_length           = rna_seq_alignment["read_length"]
  mean_fragment_length  = rna_seq_alignment["mean_fragment_length"]
  mean_dist_arg         = ""
  if paired
    files_arg = pairs.collect {|p| p[0] }.join(",") + " " + pairs.collect {|p| p[1] }.join(",")
    mean_dist           = mean_fragment_length - (read_length*2)
    mean_dist_arg       = "-r #{mean_dist_arg}"
  else
    files_arg = pairs.collect {|p| p[0] }.join(",")
  end
end
#-r = mean distance between ends of paired reads. e.g. read_length=36 & fragment_mean=300 => dist = 228
#-p = num_threads

`tophat -p #{NUM_THREADS} #{mean_dist_arg} #{files_arg} -G #{USEFUL_BED_FILES}/mm9.ucsc.genes.gtf --library-type fr-unstranded`