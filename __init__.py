import os
import sys
from .default_values import GITHUB_FOLDER
sys.path.append(GITHUB_FOLDER)
print(GITHUB_FOLDER)

import pytools.systemtools as systemtools
import pytools.filetools as filetools
import pytools.tabletools as tabletools
import pytools.timetools as timetools