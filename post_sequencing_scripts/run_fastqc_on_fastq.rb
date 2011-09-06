#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
def alive?(pid)
  pid = Integer("#{ pid }")
  begin
    Process::kill 0, pid
    true
  rescue Errno::ESRCH
    false
  end
end

require 'constants'
files = Dir.glob("#{FASTQ_CHIP_FOLDER}/**/*_fastq.*") + Dir.glob("#{FASTQ_RNA_SEQ_FOLDER}/**/*_fastq.*")
forks = []
for path in files
  forks << fork do
    name                          = File.basename(path).gsub(".gz","").gsub("_fastq.txt","")
    user                          = name.split("_")[0]
    running_file                  = running_file(name, "fastqc")
    fastqc_tmp_folder_path        = path.gsub(".txt","").gsub(".gz","") + "c"
    fastqc_tmp_zip_path           = "#{path}c.zip"
    fastqc_output_folder_path     = "#{FASTQC_FOLDER}/#{user}/#{name}_fastqc"

    next if File.exists? fastqc_output_folder_path # The file has been processed in the past
    next if File.exists? running_file #This is being processed
    puts name
    `touch #{running_file}`
    begin
      cmd = "fastqc #{path}"
      `#{cmd}`
      `mkdir -p #{FASTQC_FOLDER}/#{user}/`
      FileUtils.mv(fastqc_tmp_folder_path,fastqc_output_folder_path)
    rescue => e
      FileUtils.rm(fastqc_output_folder_path,     :force=>true)
      FileUtils.rm(fastqc_tmp_zip_path,           :force=>true)
      FileUtils.rm(fastqc_tmp_folder_path,        :force=>true)
      throw e
    ensure
      FileUtils.rm(fastqc_tmp_folder_path,        :force=>true)
      FileUtils.rm(fastqc_tmp_zip_path,           :force=>true)
      FileUtils.rm(running_file,                    :force=>true)
    end
  end
  while forks.length >= 20
    puts "#{forks.length} running. 20 max."
    Process.wait(forks[0])
    forks.delete_at(0)
  end
  
end
for fork in forks
  puts "P: Waiting for child (#{fork})"
  Process.wait(fork)
end