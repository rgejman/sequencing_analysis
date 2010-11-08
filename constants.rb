require 'rubygems'
require 'daemons'
require 'fileutils'
require 'pp'
path = File.expand_path(File.dirname(__FILE__))
CONFIG = YAML.load_file("#{path}/config.yml")


SCRIPTS_FOLDER          = CONFIG["scripts_folder"]
BASE_FOLDER             = CONFIG["base_folder"]
GENOMICS_FOLDER         = CONFIG["genomics_folder"]
RAW_FOLDER	            = "#{BASE_FOLDER}/raw"
FASTQ_FOLDER	        	= "#{BASE_FOLDER}/fastq"
FASTQC_FOLDER		        = "#{BASE_FOLDER}/fastqc"
IGVTOOLS_OUTPUT_FOLDER	= "#{BASE_FOLDER}/coverage" 
ALIGNMENTS_FOLDER       = "#{BASE_FOLDER}/alignments"
CEAS_FOLDER             = "#{BASE_FOLDER}/ceas"
LOG_FOLDER              = "#{BASE_FOLDER}/log"
TMP_FOLDER              = "#{BASE_FOLDER}/tmp"
MACS_FOLDER             = "#{BASE_FOLDER}/macs"
QUEST_FOLDER            = "#{BASE_FOLDER}/quest"
COMPOSITE_PLOTS_FOLDER  = "#{BASE_FOLDER}/composite_plots"
TOPHAT_FOLDER           = "#{BASE_FOLDER}/tophat"
MAX_FORKS               = CONFIG["max_forks"]

GENOMES_FOLDER          = "#{GENOMICS_FOLDER}/igv_tools_genomes"
CEAS_ANNOTATION_TABLES  = "#{GENOMICS_FOLDER}/ceas_annotation_tables"
QUEST_GENOME_TABLES     = "#{GENOMICS_FOLDER}/QuEST_genome_table_files"
USEFUL_BED_FILES        = "#{GENOMICS_FOLDER}/useful_bed_files"

GENOME                  = "mm9"

MYSQL_DB                = CONFIG["mysql_db"]
MYSQL_USER              = CONFIG["mysql_user"]
MYSQL_PASS              = CONFIG["mysql_pass"]
MYSQL_HOST              = CONFIG["mysql_host"]


def s_dt
  Time.now.strftime("%Y_%m_%d-%H_%M")
end

def running_file(name, tool)
  return "#{LOG_FOLDER}/#{name}.#{tool}.running"
end

def call(script_name, options)
  Daemons.call(options) do
    Dir.chdir(SCRIPTS_FOLDER)
    `ruby #{SCRIPTS_FOLDER}/#{script_name}`
  end
end