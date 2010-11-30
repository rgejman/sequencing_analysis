#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 24

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM sequencing_run WHERE run_at IS NOT NULL and illumina_run_id IS NOT NULL ORDER BY run_at desc")
res.each_hash do |sequencing_run|
  run_id                        = sequencing_run["illumina_run_id"]
  running_file                  = running_file(run_id, "BclToQseq")
  next if File.exists? running_file
  intensities_folder            = "#{RAW_FOLDER}/#{run_id}/Data/Intensities"
  base_calls_folder             = "#{intensities_folder}/BaseCalls"
  output_folder                 = base_calls_folder
  samples = {}
  samples_res = conn.query("SELECT * FROM sequencing_samples WHERE sequencing_run_id = '#{run_id}' ORDER BY lane desc")
  date    = sequencing_run["run_at"][0,10].gsub("-","_")
  nsamples_already_converted = 0
  samples_res.each_hash do |sample|
    lane = sample["lane"].to_i
    user    = sample["user"]
    name    = sample["name"]
    if user.downcase == "control"
      name += "_#{run_id}_#{lane}"
      sample["qseq_file"] = "#{QSEQ_FOLDER}/#{user}_#{date}_#{name}_qseq.txt"
    else
      sample["qseq_file"] = "#{QSEQ_FOLDER}/#{user}_#{name}_#{date}_qseq.txt"
    end
    
    nsamples_already_converted +=1 if File.exists sample["qseq_file"]
    
    samples[lane] = sample
  end
  next if nsamples_already_converted == 8 #All samples have already been converted to qseq.

  `touch #{running_file}`
  begin
    cmd = "setupBclToQseq.py -i #{base_calls_folder} -p #{intensities_folder} -o #{output_folder} --in-place --overwrite"
    puts cmd
    `#{cmd}`
    Dir.chdir(output_folder)
    cmd = "make -j #{NUM_THREADS}"
    puts cmd
    `cmd`
    for lane in sample.keys
      sample        = samples[lane]
      qseq_files    = Dir.glob("s_#{lane}_*_qseq.txt")
      qseq_file = sample["qseq_file"]
      cat_cmd = "cat #{qseq_files.join(" ")} > #{qseq_file}"
    end
  rescue => e
    throw e
  ensure
    FileUtils.rm(running_file,    :force=>true)
    conn.close
  end
  break
end