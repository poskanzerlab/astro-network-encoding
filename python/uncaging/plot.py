"""
Common plotting routines for validation figures

Designed to accompany the paper "Network-level encoding of local
neurotransmitters in cortical astrocytes" (DOI: TODO)
"""


## ------------------------------------------------------------------------ ##
## Imports
## ------------------------------------------------------------------------ ##

## ** Typing

from matplotlib.pyplot import (
    Axes,
)
from numpy.typing import (
    ArrayLike,
)


## ** Other imports

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


## ------------------------------------------------------------------------ ##
## Helper functions
## ------------------------------------------------------------------------ ##

def plot_boot( ax: Axes,
               boots: ArrayLike,
               value: float,
               color: str = 'k',
               significance: float = 0.05,
               **kwargs ):
    """
    Plot a validation of bootstrap results into the axes `ax`
    
    Parameters
    ----------
    ax : Axes
        axes to plot into
    boots : ArrayLike
        values from individual bootstrap iterations
    value : float
        actual value of the statistic
    color : str, optional
        `matplotlib` color spec for the plot (default: 'k' [black])
    significance : float, optional
        alpha cutoff for showing confidence bounds (default: 0.05)
    **kwargs : dict
        passed to `ax.hist`
    """
    
    ax.hist( boots,
             density = True,
             color = color,
             alpha = 0.5,
             **kwargs )
    
    yl = plt.ylim()
    y_std = -0.4
    y_ci = -0.8
    
    std = np.std( boots )
    low = np.quantile( boots, significance / 2. )
    high = np.quantile( boots, 1. - (significance / 2.) )
    middle = np.quantile( boots, 0.5 )
    
    ax.plot( value, y_std, f'{color}.',
             markersize = 8,
             alpha = 0.5 )
    ax.plot( [value - std, value + std], [y_std, y_std], f'{color}-',
             linewidth = 2,
             alpha = 0.5 )
    
    ax.plot( middle, y_ci, f'{color}.',
             markersize = 4,
             alpha = 0.5 )
    ax.plot( [low, high], [y_ci, y_ci], f'{color}-',
             linewidth = 1,
             alpha = 0.5 )
    
def plot_perm( ax: Axes,
               perms: ArrayLike,
               value: float,
               null: float = 0.,
               color: str = 'k',
               **kwargs ):
    """
    Plot a validation of permutation test results into an axes `ax`
    
    Parameters
    ----------
    ax : Axes
        axes to plot into
    perms : ArrayLike
        values from individual permutation iterations
    value : float
        actual value of the statistic
    null : float, optional
        the appropriate expected null value for the current setup
    color : str, optional
        `matplotlib` color spec for the plot (default: 'k' [black])
    **kwargs : dict
        passed to `ax.hist`
    """
    
    ax.hist( perms,
             density = True,
             color = color,
             alpha = 0.5,
             label = 'Permuted',
             **kwargs )
    
    yl = plt.ylim()
    
    ax.plot( [null, null], yl, f'{color}--', linewidth = 1 )
    ax.plot( [value, value], yl, f'{color}-', linewidth = 2, label = 'Observed' )
    
    ax.set_ylim( yl )