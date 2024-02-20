"""
Common statistical testing routines

Designed to accompany the paper "Network-level encoding of local
neurotransmitters in cortical astrocytes" (DOI: TODO)
"""


## ------------------------------------------------------------------------ ##
## Imports
## ------------------------------------------------------------------------ ##

## ** Typing

from typing import (
    Protocol,
)
from pandas import (
    DataFrame,
)
from numpy.typing import (
    ArrayLike,
    NDArray,
)

# Protocols (for typing)

class DFStatistic( Protocol ):
    def __call__( self,
                  df: DataFrame ) -> float:
        return super().__call__( df )


## ** Other imports

import numpy as np
import pandas as pd

from tqdm import tqdm

## ** Local imports

from . import events


## ------------------------------------------------------------------------ ##
## Helper functions
## ------------------------------------------------------------------------ ##

def perm_p(
        n0: int,
        n: int
    ) -> float:
    """
    Compute conservative permutation test p-value

    The bias introduced by adding one to the numerator and denominator prevents
    p-values of 0, which can lead to downstream issues with Bayesian reasoning,
    multiple test corrections, etc.

    For details, see TODO citation
    
    Parameters
    ----------
    n0 : int
        number of rejections under the permuted null
    n : int
        total number of permutations
    
    Returns
    -------
    float
        the (appropriately) biased permutation p-value
    """
    return (n0 + 1) / (n + 1)

def bootstrap_sample_hierarchical(
        df: DataFrame,
        hierarchy: 'list[str]',
        fix: 'list[str]' = []
    ) -> DataFrame:
    """
    Perform sampling from `df` with resampling at each level in `hierarchy`
    
    Parameters
    ----------
    df : DataFrame
        the input data
    hierarchy : list[str]
        ordering of the levels of the resampling hierarchy, most coarse first to
        most fine last
    fix : list[str], optional
        list of levels to not resample (default: []); this makes it so that the
        *hierarchy* at this level is maintained, but there is no *randomness* in
        the sampling from groupings at this level
    
    Returns
    -------
    DataFrame
        the bootstrapped sample, in the same format as `df`
    """

    if len( hierarchy ) == 0:
        # Base case:
        # Nothing left to bootstrap over
        return df
    
    level = hierarchy[0]

    # Determine which values at the current level to sample
    keys = df[level].unique()
    if level in fix:
        # Level is fixed; use everything as-is
        selected = keys
    else:
        # Level is not fixed; sample with replacement
        selected = np.random.choice( keys, len( keys ) )
    
    # Recursion:
    # The sample in *each value* at the current level is formed by bootstrapping
    #   sampling over all of the data for that value with the finer part of the
    #   hierarchy
    return pd.concat( [ bootstrap_sample_hierarchical( df[df[level] == x],
                                                       hierarchy[1:],
                                                       fix = fix )
                        for x in selected ] )

def bootstrap_stat_hierarchical(
        f: DFStatistic,
        df: DataFrame,
        hierarchy: 'list[str]',
        n: int = 200,
        verbose: bool = False,
        **kwargs
    ) -> NDArray:
    """
    Perform hierarchical bootstrapping of the statistic `f`
    
    Parameters
    ----------
    f : DFStatistic
        the statistic: a function that takes in a `DataFrame` (bootstrap sample)
        and returns a scalar
    df : DataFrame
        the original data
    hierarchy : list[str]
        a list specifying the hierarchical structure over which to bootstrap;
        the first entry is the top level of the hierarchy
    n : int, optional
        the number of bootstrap samples to perform (default: 200)
    verbose : bool, optional
        if `True`, show a progress bar
    **kwargs : dict, optional
        passed to `bootstrap_sample_hierarchical
    
    Returns
    -------
    NDArray
        the values of `f` for each bootstrap iteration
    """
    
    ret = np.zeros( (n,) )
    
    if verbose:
        it = tqdm( range( n ) )
    else:
        it = range( n )
        
    for i in it:
        df_boot = bootstrap_sample_hierarchical( df, hierarchy, **kwargs )
        ret[i] = f( df_boot )
        
    return ret

def perm_diff_stat_df(
        f: DFStatistic,
        df: DataFrame,
        key: str,
        n: int = 1000,
        verbose: bool = False
    ) -> NDArray:
    """
    Perform permutation testing of `f` by shuffling `key` within the data `df`
    
    Parameters
    ----------
    f : DFStatistic
        the statistic, a function that takes in a `DataFrame` (permuted sample)
        and returns a scalar
    df : DataFrame
        the original data
    key : str
        the key within `df` to permute each iteration
    n : int, optional
        the number of permuted samples to compute (default: 1000)
    verbose : bool, optional
        if `True`, show a progress bar

    Returns
    -------
    NDArray
        the values of `f` for each permutation iteration
    """
    
    ret = np.zeros( (n,) )
    
    if verbose:
        it = tqdm( range( n ) )
    else:
        it = range( n )
        
    for i in it:
        df_perm = df.copy()
        df_perm[key] = np.random.permutation( df_perm[key] )
        ret[i] = f( df_perm )
    
    return ret

def perm_stat_ts(
        f: DFStatistic,
        df: DataFrame,
        window: events.Window,
        t_key: str = 'start_time_rel',
        shift_level: str = 'cell_global_all',
        n: int = 200
    ) -> NDArray:
    """
    TODO

    Parameters
    ----------
    f : DFStatistic
        the statistic, a function that takes in a `DataFrame` (permuted sample)
        and returns a scalar
    df : DataFrame
        the original data
    window : events.Window
        the window within which to randomly shift the data; see
        `events.perm_events_shift_group`
    t_key : str
        the key in `df` that corresponds to time (default: 'start_time_rel');
        see `events.perm_events_shift_group`
    shift_level : str
        the key in `df` to group events together when shifting as a block
        (default: 'cell_global_all'); see `events.perm_events_shift_group`
    n : int
        the number of permutation iterations (default: 200)
    
    Returns
    -------
    NDArray
        the values of `f` for each shift-permutation iteration
    """
    
    ret = np.zeros( (n,) )
    
    for i in range( n ):
        df_perm = events._perm_events_shift_group(
            df,
            window,
            t_key = t_key,
            level = shift_level
        )
        ret[i] = f( df_perm )                           
    
    return ret