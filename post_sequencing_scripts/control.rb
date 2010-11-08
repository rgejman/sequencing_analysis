#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'

options = {
    :ARGV       => ['start'], #, '-f', '--', 'param_for_myscript'
    :dir_mode   => :normal,
    :dir        => LOG_FOLDER,
    :multiple   => true,
    :ontop      => false,
    :mode       => :load,
    :backtrace  => false,
    :monitor    => false,
    :log_output => true
  }

programs = ["run_fastqc.rb", "run_alignment.rb", "run_sort_sam.rb",
              "run_bam_index.rb", "run_igvtools.rb", "run_make_bam.rb",
                "run_macs_single.rb","run_macs_pair.rb", "run_quest.rb",
                "run_tophat.rb"] # #,"run_ceas.rb" => Need to make sure is integrated.
              #,"run_make_bed.rb","run_make_wig.rb" # unclear if these are necessary
              
def forks
  Dir.entries("#{LOG_FOLDER}/").select{|e| e =~ /\.pid/ }.length
end

loop do
  for p in programs
    while forks >= MAX_FORKS
      puts "#{forks} running. #{MAX_FORKS} max."
      sleep(30)
    end
    puts "#{forks} running. #{MAX_FORKS} max."
    call(p, options)
    puts "Called #{p}"
  end
  sleep(30)
end