import os

PIPELINE_FOLDER = "/home/upmc/Documents/Variant_Discovery_Pipeline/"
default_config_filename = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_project_options.txt")
#default_sample_filename = os.path.join(PIPELINE_FOLDER, "debug_sample_list.tsv")
default_sample_filename = os.path.join(PIPELINE_FOLDER, "debug_sample_list.tsv")
pipeline_dna_callers = ['missing']#['varscan', 'somaticsniper', 'muse', 'mutect', 'strelka']
pipeline_rna_callers = []
pipeline_copynumber_callers = []
somaticseq_callers = ['table']