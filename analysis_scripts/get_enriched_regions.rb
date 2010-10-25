$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
ARGS_REQUIRED = 2
unless ARGV.length == ARGS_REQUIRED
  puts "Usage:\t\tget_enriched_regions.rb region_file enrichment_file"
  files = Dir.glob("/Users/rgejman/genomics/useful_bed_files/*.bed").collect {|f| f.split("/").last}
  puts ""
  puts "Region files:\t#{files.join("\n\t\t")}"
  exit
end

region          = ARGV[0]
enrichment_file = ARGV[1]

puts `intersectBed -u -wa -a #{region_file} -b #{enrichment_file}` # write each "A" line only once, no matter how many matches.