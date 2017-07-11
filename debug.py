

PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
default_config_filename = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")

default_sample_filename = os.path.join(PIPELINE_DIRECTORY, "rna_sample_list.tsv")

pipeline_dna_callers = []
pipeline_rna_callers = ['haplotypecaller']
pipeline_copynumber_callers = []
somaticseq_callers = ['table'] # Different from the others. should be ['mode', 'truthset path', 'classifier path']