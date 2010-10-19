#!/usr/bin/env ruby -KU
require 'constants'

options = {
    :ARGV       => ['start'], #, '-f', '--', 'param_for_myscript'
    :dir_mode   => :normal,
    :dir        => LOG_FOLDER,
    :multiple   => true,
    :ontop      => false,
    :mode       => :load,
    :backtrace  => true,
    :monitor    => false,
    :log_output => true
  }

forks = []
programs = ["run_fastqc.rb", "run_alignment.rb", "run_sort_sam.rb",
              "run_bam_index.rb", "run_igvtools.rb", "run_make_bam.rb", "run_macs.rb", "run_quest.rb"] # #,"run_ceas.rb" => Need to make sure is integrated.
              #,"run_make_bed.rb","run_make_wig.rb" # unclear if these are necessary

loop do
  for p in programs
    while forks.length >= MAX_FORKS
      sleep(5)
      forks.delete_if {|f| !f.running? }
    end
    forks << call(p, options)
  end
end