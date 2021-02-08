#Module containing the imported libraries

import random
import time
import csv
#from itertools import izip
import itertools
import numpy as np
from collections import Counter
from sklearn.preprocessing import scale
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ExpSineSquared
from scipy.stats import norm
from scipy.optimize import fsolve
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
