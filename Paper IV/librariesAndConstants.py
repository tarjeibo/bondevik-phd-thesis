#Libraries
import csv
import os
import numpy as np
from cmath import sqrt, log10
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from decimal import Decimal

#Physical constants
e = 1.60217662e-19
epsilon_0 = 8.85e-12
j = sqrt(-1)
k_B = 1.38064852e-23

#BaZrO3 material specific constants
epsilon_r = 75
a_0 = 4.235e-10
c_acc = 0.05
