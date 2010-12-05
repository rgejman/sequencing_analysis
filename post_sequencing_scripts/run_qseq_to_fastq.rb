#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
samples_res = conn.query("SELECT * FROM sequencing_samples,sequencing_run where sequencing_run_id = sequencing_run.id and user != 'Control' and user != 'control'")

samples_res.each_hash do |sample|
  date          = sample["run_at"][0,10].gsub("-","_")
  lane          = sample["lane"].to_i
  user          = sample["user"].capitalize
  name          = sample["name"]
  type          = sample["type"]
  run_id        = sample["sequencing_run_id"]
  base_file     = sample_filebase(run_id, date, lane, user, name) + "_fastq.txt"
  qseq_file     = "#{QSEQ_FOLDER}/#{user}/" + base_file
  folder        = (type == "chip" ? "#{FASTQ_CHIP_FOLDER}" : "#{FASTQ_RNA_SEQ_FOLDER}") + "/#{user}"
  fastq_file    = "#{folder}/#{base_file}"
  tmp_file      = "#{TMP_FOLDER}/#{base_file}"
  running_file  = running_file(base_file, "qseq_to_fastq")
  puts fastq_file
  puts qseq_file
  next unless File.exists? qseq_file
  next if File.exists? fastq_file
  next if File.exists? running_file
  `touch #{running_file}`  
  
  Dir.chdir(QSEQ_FOLDER)
  begin
    cat  = "cat #{file}"
    awk1 = "gawk '{gsub(/\\./, \"N\", $9);print}'"
    awk2 = "awk '{print \"@\"$1\":#{base}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\"; print $9; print \"+\"$1\":#{base}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\" ; print $10}'"
    out  = "> #{tmp_file}"
    `#{cat} | #{awk1} | #{awk2} #{out}`
    `mkdir -p #{folder}`
    FileUtils.mv(tmp_file, fastq_file)
  rescue => e
    FileUtils.rm(fastq_file,     :force=>true)
    FileUtils.rm(tmp_file,       :force=>true)
    throw e
  ensure
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
  
end