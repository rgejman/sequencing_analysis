#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 24

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM sequencing_run WHERE run_at IS NOT NULL and illumina_run_id IS NOT NULL ORDER BY run_at desc")
res.each_hash do |sequencing_run|
  seq_run_id                    = sequencing_run["id"]
  run_id                        = sequencing_run["illumina_run_id"]
  paired                        = sequencing_run["paired"].to_i == 1
  running_file                  = running_file(run_id, "BclToQseq")
  next if File.exists? running_file
  intensities_folder            = "#{RAW_FOLDER}/#{run_id}/Data/Intensities"
  base_calls_folder             = "#{intensities_folder}/BaseCalls"
  samples = {}
  samples_res = conn.query("SELECT * FROM sequencing_samples WHERE sequencing_run_id = '#{seq_run_id}' ORDER BY lane desc")
  date    = sequencing_run["run_at"][0,10].gsub("-","_")
  nsamples_not_converted = 0
  samples_res.each_hash do |sample|
    lane        = sample["lane"].to_i
    user        = sample["user"].capitalize
    name        = sample["name"]
    type        = sample["type"]
    base        = sample_filebase(run_id, date, lane, user, name)
    next if sample["post_process"].to_i == 0

    qseq_base     = "#{QSEQ_FOLDER}/#{user}/"
    chip_base     = "#{FASTQ_CHIP_FOLDER}/#{user}/"
    rna_seq_base  = "#{FASTQ_RNA_SEQ_FOLDER}/#{user}/"
    hic_base      = "#{FASTQ_HIC_FOLDER}/#{user}/"

    if paired
      sample["qseq_file_1"]   = base + "_1_qseq.txt"
      sample["qseq_file_2"]   = base + "_2_qseq.txt"
      fq1 = sample["qseq_file_1"] + "_1_fastq.txt"
      fq2 = sample["qseq_file_2"] + "_2_fastq.txt"
      next if File.exists?("#{qseq_base}" + sample["qseq_file_1"]) and File.exists?("#{qseq_base}" + sample["qseq_file_2"])
      next if File.exists?("#{chip_base}" + fq1) and File.exists?("#{chip_base}" + fq2)
      next if File.exists?("#{rna_seq_base}" + fq1) and File.exists?("#{rna_seq_base}" + fq2)
      next if File.exists?("#{hic_base}" + fq1) and File.exists?("#{hic_base}" + fq2)
      nsamples_not_converted += 1 if user.downcase != "control"
    else
      sample["qseq_file"]     = base + "_qseq.txt"
      fq = base + "_fastq.txt"
      q = qseq_base + sample["qseq_file"]
      c = chip_base + fq
      r = rna_seq_base + fq
      c_gz = "#{FASTQ_CHIP_FOLDER}/#{user}/#{fq}.gz"
      r_gz = "#{FASTQ_RNA_SEQ_FOLDER}/#{user}/#{fq}.gz"
      next if File.exists? q or File.exists? r or File.exists? c or File.exists? c_gz or File.exists? r_gz
      nsamples_not_converted += 1 if user.downcase != "control"
    end
    samples[lane] = sample
  end
  next if nsamples_not_converted == 0

  `touch #{running_file}`
  begin
    cmd = "setupBclToQseq.py -b #{base_calls_folder} --in-place --overwrite"
    puts cmd
    `#{cmd}`
    Dir.chdir(base_calls_folder)
    cmd = "make -j #{NUM_THREADS}"
    puts cmd
    `#{cmd}`
    forks = []
    for lane in samples.keys
      sample        = samples[lane]
      ##### WE DO NOT CONVERT CONTROL LANES TO QSEQ. WASTE OF SPACE! 
      next if sample["user"].downcase == "control"

      puts "P: Cat'ing lane #{lane}"
      if paired
        qseq_files_1    = Dir.glob("s_#{lane}_1_*_qseq.txt").sort
        qseq_files_2    = Dir.glob("s_#{lane}_2_*_qseq.txt").sort
        qseq_file_1     = sample["qseq_file_1"]
        qseq_file_2     = sample["qseq_file_2"]
        
        qseq_filepath_1 = "#{QSEQ_FOLDER}/#{sample['user'].capitalize}/" + qseq_file_1
        qseq_filepath_2 = "#{QSEQ_FOLDER}/#{sample['user'].capitalize}/" + qseq_file_2
        
        tmp_filepath_1  = "#{TMP_FOLDER}/" + qseq_file_1
        tmp_filepath_2  = "#{TMP_FOLDER}/" + qseq_file_2
        
        ## Concatenate all the tiles and strip out the "failed" reads (to save on space and aligning later, etc)
        cat_cmd_1 = "cat #{qseq_files_1.join(" ")} > #{tmp_filepath_1}"
        cat_cmd_2 = "cat #{qseq_files_2.join(" ")} > #{tmp_filepath_2}"
        puts "C: #{cat_cmd_1}"
        puts "C: #{cat_cmd_2}"
        forks << fork do
          `#{cat_cmd_1}`
          `#{cat_cmd_2}`
          `mkdir -p #{QSEQ_FOLDER}/#{sample['user'].capitalize}/`
          FileUtils.mv(tmp_filepath_1, qseq_filepath_1)
          FileUtils.mv(tmp_filepath_2, qseq_filepath_2)
        end
      else
        qseq_files    = Dir.glob("s_#{lane}_*_qseq.txt")
        qseq_file     = sample["qseq_file"]
        qseq_filepath = "#{QSEQ_FOLDER}/#{sample['user'].capitalize}/" + sample["qseq_file"]
        tmp_filepath  = "#{TMP_FOLDER}/" + sample["qseq_file"]
        ## Concatenate all the tiles and strip out the "failed" reads (to save on space and aligning later, etc)
        cat_cmd = "cat #{qseq_files.join(" ")} | grep -e \"1$\" > #{tmp_filepath}"
        puts "C: #{cat_cmd}"
        forks << fork do
          `#{cat_cmd}`
          `mkdir -p #{QSEQ_FOLDER}/#{sample['user'].capitalize}/`
          FileUtils.mv(tmp_filepath, qseq_filepath)
        end

      end
    end
    for fork in forks
      puts "P: Waiting for child (#{fork})"
      Process.wait(fork)
    end
  rescue => e
    throw e
  ensure
    FileUtils.rm(running_file,    :force=>true)
    conn.close
  end
  break
end