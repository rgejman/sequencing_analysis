#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'open3'
EXTEND_LENGTH = 100
WINDOW_SIZE   = 25
files = Dir.glob("#{ALIGNMENTS_FOLDER}/**/*.sorted.bam")

for file in files
  pid = fork do
    base                  = file.split("/").last.gsub('.bam','') #.bam is added automatically
    user                  = base.split("_").first
    running_file          = running_file(base, "coverage")
    tmp_file1             = "#{TMP_FOLDER}/#{base}.cov.tdf"
    tmp_file2             = "#{TMP_FOLDER}/#{base}.wig"
    tmp_file3             = tmp_file2.gsub(".sorted",".normalized.sorted")
    output_file_1         = "#{COVERAGE_FOLDER}/#{user}/#{base}.cov.tdf"
    output_file_2         = "#{WIG_FOLDER}/#{user}/#{base}.#{WINDOW_SIZE}.wig"
    input_path            = file
    genome = "mm9"
    if file =~ /A549/ 
      genome = "hg19"
    end
    next if File.size(file) < 1024*1024 # 1 megabyte
    next if File.exists?(output_file_1) and File.exists?(output_file_2) #Already processed
    next if File.exists? running_file #Being processed
    puts file
    `touch #{running_file}`
    begin
      `igvtools count -e #{EXTEND_LENGTH} -w #{WINDOW_SIZE} "#{input_path}" "#{tmp_file1},#{tmp_file2}" "#{GENOMES_FOLDER}/#{GENOME}.genome"`
      num_mapped_reads = 0
      c =  "bamstats #{file}"
      puts c
      Open3.popen3(c) { |stdin, stdout, stderr|
        while(line = stdout.gets)
          rd = line if line =~ /Mapped/
        end
        num_mapped_reads = rd.split(/\s+/)[3].to_i
      }
      norm = 1000000.0/num_mapped_reads.to_f
      puts "#{file}: #{num_mapped_reads} reads. Will normalize by #{norm}."
      `normalizeWig #{norm} #{tmp_file2}`
      `mkdir -p #{COVERAGE_FOLDER}/#{user}`
      `mkdir -p #{WIG_FOLDER}/#{user}`
      FileUtils.mv(tmp_file1, output_file_1)
      FileUtils.mv(tmp_file3, output_file_2)
    rescue => e
      FileUtils.rm(output_file_1,     :force=>true)
      FileUtils.rm(output_file_2,     :force=>true)
      throw e
    ensure
      FileUtils.rm(tmp_file1,        :force=>true)
      FileUtils.rm(tmp_file2,        :force=>true)
      FileUtils.rm(tmp_file3,        :force=>true)
      FileUtils.rm(running_file,    :force=>true)
    end
  end
end