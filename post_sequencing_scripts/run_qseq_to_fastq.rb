#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
samples_res = conn.query("SELECT * FROM sequencing_samples,sequencing_run where sequencing_run_id = sequencing_run.id and user != 'Control' and user != 'control'")
forks = []
samples_res.each_hash do |sample|
  date          = sample["run_at"][0,10].gsub("-","_")
  lane          = sample["lane"].to_i
  user          = sample["user"].capitalize
  name          = sample["name"]
  type          = sample["type"].downcase
  run_id        = sample["sequencing_run_id"]
  paired        = sample["paired"].to_i == 1
  base_file     = sample_filebase(run_id, date, lane, user, name)
  if paired
    qseq_1_file     = "#{QSEQ_FOLDER}/#{user}/" + base_file + "_1_qseq.txt"
    qseq_2_file     = "#{QSEQ_FOLDER}/#{user}/" + base_file + "_2_qseq.txt"
  else
    qseq_file     = "#{QSEQ_FOLDER}/#{user}/" + base_file + "_qseq.txt"
  end
  folder = case type.downcase
  when "chip"
    FASTQ_CHIP_FOLDER
  when "rna_seq"
    FASTQ_RNA_SEQ_FOLDER
  when "hic"
    FASTQ_HIC_FOLDER
  end
  running_file  = running_file(base_file, "qseq_to_fastq")
  folder        += "/#{user}"
  if paired
    fastq_1_file    = "#{folder}/#{base_file}_1_fastq.txt"
    fastq_2_file    = "#{folder}/#{base_file}_2_fastq.txt"
    tmp_1_file      = "#{TMP_FOLDER}/#{base_file}_1_fastq.txt"
    tmp_2_file      = "#{TMP_FOLDER}/#{base_file}_2_fastq.txt"
    next unless File.exists? qseq_1_file and File.exists? qseq_2_file
    next if File.exists? fastq_1_file and File.exists fastq_2_file
  else
    fastq_file    = "#{folder}/#{base_file}_fastq.txt"
    tmp_file      = "#{TMP_FOLDER}/#{base_file}_fastq.txt"
    next unless File.exists? qseq_file
    next if File.exists? fastq_file
  end
  next if File.exists? running_file
  next if sample["post_process"].to_i == 0
  `touch #{running_file}`
  forks << fork do
    if paired
      puts qseq_1_file + " and " + qseq_2_file
    else
      puts qseq_file
    end
    Dir.chdir(QSEQ_FOLDER)
    begin
      if paired
        cat1  = "cat #{qseq_1_file}"
        cat2  = "cat #{qseq_2_file}"
        awk1 = "gawk '{gsub(/\\./, \"N\", $9);print}'"
        awk2 = "awk '{print \"@\"$1\":#{base_file}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\"; print $9; print \"+\"$1\":#{base_file}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\" ; print $10}'"
        out1  = "> #{tmp_1_file}"
        out2  = "> #{tmp_2_file}"
        
        `#{cat1} | #{awk1} | #{awk2} #{out1}`
        `#{cat2} | #{awk1} | #{awk2} #{out2}`
        
        `mkdir -p #{folder}`
        FileUtils.mv(tmp_1_file, fastq_1_file)
        FileUtils.mv(tmp_2_file, fastq_2_file)
      else
        cat  = "cat #{qseq_file}"
        awk1 = "gawk '{gsub(/\\./, \"N\", $9);print}'"
        awk2 = "awk '{print \"@\"$1\":#{base_file}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\"; print $9; print \"+\"$1\":#{base_file}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\" ; print $10}'"
        out  = "> #{tmp_file}"
        `#{cat} | #{awk1} | #{awk2} #{out}`
        `mkdir -p #{folder}`
        FileUtils.mv(tmp_file, fastq_file)
      end
    rescue => e
      if paired
        FileUtils.rm(fastq_1_file,     :force=>true)
        FileUtils.rm(fastq_2_file,     :force=>true)
        FileUtils.rm(tmp_1_file,       :force=>true)
        FileUtils.rm(tmp_2_file,       :force=>true)
      else
        FileUtils.rm(fastq_file,     :force=>true)
        FileUtils.rm(tmp_file,       :force=>true)
      end
      throw e
    ensure
      FileUtils.rm(running_file,    :force=>true)
    end  
  end
end