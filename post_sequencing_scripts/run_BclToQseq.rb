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
  running_file                  = running_file(run_id, "BclToQseq")
  next if File.exists? running_file
  intensities_folder            = "#{RAW_FOLDER}/#{run_id}/Data/Intensities"
  base_calls_folder             = "#{intensities_folder}/BaseCalls"
  output_folder                 = base_calls_folder
  samples = {}
  samples_res = conn.query("SELECT * FROM sequencing_samples WHERE sequencing_run_id = '#{seq_run_id}' ORDER BY lane desc")
  date    = sequencing_run["run_at"][0,10].gsub("-","_")
  nsamples_not_converted = 0
  samples_res.each_hash do |sample|
    lane        = sample["lane"].to_i
    user        = sample["user"].capitalize
    name        = sample["name"]
    base        = sample_filebase(run_id, date, lane, user, name)
    sample["qseq_file"]   = base + "_qseq.txt"
    fq  = base + "_fastq.txt"
    q = "#{QSEQ_FOLDER}/#{user}/" + sample["qseq_file"]
    c = "#{FASTQ_CHIP_FOLDER}/#{user}/#{fq}"
    r = "#{FASTQ_RNA_SEQ_FOLDER}/#{user}/#{fq}"
    c_gz = "#{FASTQ_CHIP_FOLDER}/#{user}/#{fq}.gz"
    r_gz = "#{FASTQ_CHIP_FOLDER}/#{user}/#{fq}.gz"
    
    if File.exists? q or File.exists? r or File.exists? c or File.exists? c_gz or File.exists? r_gz
      next
    else
      nsamples_not_converted += 1 if user.downcase != "control"
    end
    samples[lane] = sample
  end
  next if nsamples_not_converted == 0

  `touch #{running_file}`
  begin
    cmd = "setupBclToQseq.py -b #{base_calls_folder} -p #{intensities_folder} -o #{output_folder} --in-place --overwrite"
    puts cmd
    `#{cmd}`
    Dir.chdir(output_folder)
    cmd = "make -j #{NUM_THREADS}> /dev/null 2>&1"
    puts cmd
    `#{cmd}`
    forks = []
    for lane in samples.keys
      sample        = samples[lane]
      ##### WE DO NOT CONVERT CONTROL LANES TO QSEQ. WASTE OF SPACE! 
      next if sample["user"].downcase == "control"
      
      puts "P: Cat'ing lane #{lane}"
      
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