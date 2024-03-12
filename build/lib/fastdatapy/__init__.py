#Import the used libraries

import re
import urllib
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
import os
import datetime
import time
from array import array
import matplotlib.pyplot as plt
from pylab import *
import seaborn as sns
from matplotlib.pyplot import MultipleLocator
import scipy.interpolate as spi
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit
import statsmodels.api as sm
import shutil
from pybaselines import Baseline
from pybaselines import utils
from pybaselines.utils import gaussian
from scipy.signal import savgol_filter as sf
from decimal import *
# from matplotlib import cm
# from numba import njit, prange
from scipy.sparse.linalg import spsolve
from scipy import sparse

# from matplotlib.patches import Rectangle

lowess = sm.nonparametric.lowess




dates = np.array(['20210425', '20210426', '20210427', '20210430', '20210501', '20210502',
 '20210503', '20210504', '20210505', '20210506', '20210507', '20210508',
 '20210526'])

__all__=["fastoperations"] # 列表可以根据要导入的模块数而进行新增，列表元素是之前新建的py文件的名字

__version__ = "1.0"


