import os

GITHUB_FOLDER = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.append(GITHUB_FOLDER)
print(GITHUB_FOLDER)
import pytools.systemtools as systemtools
import pytools.filetools as filetools
import pytools.tabletools as tabletools
import pytools.timetools as timetools
import datetime
now = datetime.datetime.now
# The parent folder the pipeline will be run in.
PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
# File to save the console output to. Only used when the console output is supressed.
CONSOLE_LOG_FILE = ""
# File containing a test sample.
SAMPLE_LOG_FILE = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "sample_logV2.tsv")
README_FILE     = os.path.join(PIPELINE_DIRECTORY, "0_readme_files", "readme.{0}.txt".format(now().isoformat()))
# Whether to use backwards-compatible filenames
BACKWARDS_COMPATIBLE = False
# Whether to overwrite any existing files
FORCE_OVERWRITE = False

DEBUG = True
VERBOSE_LEVEL = 4

# ----------------------------------------------------------------------------------------------------
# ------------------------------------ Set Up Global Functions ---------------------------------------
# ----------------------------------------------------------------------------------------------------

def generateTimestamp():
	timestamp = now().isoformat()
	timestamp = timestamp.split('.')[0]
	return timestamp


def getsize(path, total = True):
	sizes = []
	if os.path.isfile(path):
		sizes += [os.path.getsize(path)]
	elif os.path.isdir(path):
		dirsize = list()
		for fn in os.listdir(path):
			new_path = os.path.join(path, fn)
			dirsize += getsize(new_path, total = False)
		sizes += dirsize
	else:
		pass
	if total: sizes = sum(sizes)
	return sizes
