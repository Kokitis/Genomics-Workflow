import os
import sys

github_folder = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(github_folder)

import pytools.tabletools as tabletools
import pytools.filetools as filetools
import pytools.timetools as timetools
import pytools.systemtools as systemtools

import gdc_api
API = gdc_api.GDCAPI()