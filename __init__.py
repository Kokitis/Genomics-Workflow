import os
import sys
GITHUB_FOLDER = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(GITHUB_FOLDER)
print(GITHUB_FOLDER)

import pytools.systemtools as systemtools
import pytools.filetools as filetools
import pytools.tabletools as tabletools
import pytools.timetools as timetools