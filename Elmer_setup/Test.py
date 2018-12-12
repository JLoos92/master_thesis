import matplotlib
#matplotlib.use('nbagg')
from pylab import *
import numpy as np
from pandas import DataFrame, Series
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from numpy import genfromtxt
import plotly as py
py.tools.set_credentials_file(username='Julius92Loos', api_key='ZdBblhNTi4fwVZUtaGwy')
from plotly import __version__
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go

# more
from scipy.stats import spearmanr
from scipy.stats import norm 
from scipy.stats import gamma
from scipy.stats import beta
from scipy.spatial import distance_matrix
from scipy import stats as sst
from scipy.stats import multivariate_normal
from scipy.spatial.distance import cdist
from IPython.display import Image
from astropy.io import ascii
from IPython.display import display
from IPython.display import SVG
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization! 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



import plotly.plotly as py
import plotly.graph_objs as go

import pandas as pd

x = np.random.randn(2000)
y = np.random.randn(2000)
iplot([go.Histogram2dContour(x=x, y=y, contours=dict(coloring='heatmap')),
       go.Scatter(x=x, y=y, mode='markers', marker=dict(color='white', size=3, opacity=0.3))], show_link=False)