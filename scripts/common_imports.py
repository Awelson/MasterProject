# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:59:15 2024

@author: aaron
"""

# Imports

import polars as pl
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
import plotnine as p9
import matplotlib.pyplot as plt
import statsmodels.api as sm
import itertools
import seaborn as sns

def rgb_to_hex(rgb_tuple):
    r, g, b = [int(c * 255) for c in rgb_tuple]
    return "#{:02x}{:02x}{:02x}".format(r, g, b)
