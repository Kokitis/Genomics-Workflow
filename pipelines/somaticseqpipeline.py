from .basepipeline import BasePipeline 
from somaticseq.somaticseq import SomaticSeq

class SomaticSeqPipeline(BasePipeline):
	def runWorkflow(self, sample, pipeline_options, pipeline_callers):
		"""
			Parameters
			----------
				sample:
				options: 
				callers: list<>
					The parameters to run somaticseq with.
					Positional Arguments:
						*mode: {'trainer', 'prediction', 'table'}
							The mode to run somaticseq in.
						*truthset: Path to a truthset or a truthset object
						*classifier: path to a somaticseq classifier
		"""

		somaticseq_mode = pipeline_callers[0]
		
		if len(pipeline_callers) == 3 or somaticseq_mode == 'trainer':
			somaticseq_truthset = pipeline_callers[1]
		else:
			somaticseq_truthset = None

		if len(pipeline_callers) == 3 or somaticseq_mode == 'prediction':
			somaticseq_classifier = pipeline_callers[-1]
		else:
			somaticseq_classifier = None

		result = SomaticSeq(
			sample = sample, 
			options = pipeline_options,
			mode = somaticseq_mode,
			truthset = somaticseq_truthset,
			classifier = somaticseq_classifier
		)
		return result

