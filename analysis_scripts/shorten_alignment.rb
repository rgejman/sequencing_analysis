#!/usr/bin/env ruby -wKU
require 'rubygems'
require 'open3'


## This script will shorten alignment files to an approximately equal # of mapped reads.
## To ensure that the shortening is done randomly w/o reading the entire alignment into memory,
## each read is challenged independently until the end of the file, yielding an approximately
## equal # of reads.
##
## Unmapped reads are discarded.
##
## This script is useful for when you have more reads in one sample than another
## and you want to compare  alignments with the same # of reads and/or to get comparable
## RPKM values.

num_alignments  = ARGV.length
alignments      = ARGV
lengths         = []
min_length      = nil

EQUIVALENCE_THRESHOLD = 0.05

for i in (0...num_alignments)
  alignment_file  = alignments[i]
  raise "File could not be found: #{alignment_file}" unless File.exists? alignment_file
  l = `samtools flagstat #{alignment_file}`.split("\n")[3].split(/\s/)[0].to_i
  lengths << l
  min_length = l if min_length.nil? or l < min_length
end
puts "*** Shrinking all files to #{min_length} ***"

for i in (0...num_alignments)
  alignment_file  = alignments[i]
  l               = lengths[i]
  if l >= (1-EQUIVALENCE_THRESHOLD) * min_length and l <= (1 + EQUIVALENCE_THRESHOLD) * min_length
    # This alignment is within 5% of the minimum size.
    puts "#{alignment_file} #{lengths[i].to_f / min_length.to_f*100}% of min size."
  elsif l > min_length
    p = min_length.to_f / l.to_f
    output_file_tokens = alignment_file.split(".")
    output_file = output_file_tokens.shift + ".approx_#{min_length}." + output_file_tokens.join(".").gsub(".bam","")
    Open3.popen3("bamtools filter -in #{alignment_file} -isMapped true | samtools view -h -") do |i_stdin, i_stdout,i_stderr|
      Open3.popen3("samtools view -hbS - | samtools sort - #{output_file}") do |o_stdin, o_stdout,o_stderr|
        t = Thread.new(o_stdout, o_stderr) do |o_o,o_e|
           while (line = o_o.gets)
             puts line
           end
           while (line = o_e.gets)
             puts line
           end
         end
        while (line = i_stdout.gets)
          if line[0,1] != "@"
            next if rand() > p
          end
          o_stdin.print line
        end
        o_stdin.close
        t.join
      end
    end
    `bamtools index -in #{output_file}.bam`
  else
    raise "Error: #{alignment_file} cannot be smaller than min_length. #{l} #{min_length}"
  end
end