"""
Miscellaneous utilities to support uncaging analyses

Designed to accompany the paper "Network-level encoding of local
neurotransmitters in cortical astrocytes" (DOI: TODO)
"""


## ------------------------------------------------------------------------ ##
## Imports
## ------------------------------------------------------------------------ ##

from typing import (
    Union,
)
from pandas import (
    DataFrame,
)
from numpy.typing import (
    ArrayLike,
    NDArray,
)

import yaml

from pathlib import Path

import numpy as np
import pandas as pd


## ------------------------------------------------------------------------ ##
## Helper functions
## ------------------------------------------------------------------------ ##

def _merge_dicts( *dicts: dict ) -> dict:
    """
    Returns a new `dict` by sequentially overwriting with each input arg
    """
    ret = dict()
    for d in dicts:
        for k, v in d.items():
            ret[k] = v
    return ret

def _zip_by_dataset(
        headers, # TODO Type annotations
        events,
    ):
    """
    TODO
    """
    ret = []
    for dataset_id, dataset_events in events.groupby( 'dataset_id' ):
        dataset_header = headers.get( dataset_id, None )
        ret.append( (dataset_header, (dataset_id, dataset_events)) )
    return ret

def p_text( p: float ) -> str:
    """
    Figure significance-shorthand string for p-value `p`
    """
    if p < 0.0001:
        return '#'
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    else:
        return 'n.s.'

def best_parallel_workers( max_workers: int,
                           n_total: int ) -> 'tuple[int, int]':
    """
    Find the right number of parallel workers to use

    Given that there are `n_total` jobs to perform and `max_workers` available
    workers, there can arise situations where using all of the available workers
    is not optimal because of integer division arithmetic. This streamlines this
    computation.

    Parameters
    ----------
    max_workers : int
        the maximum number of available parallel workers
    n_total : int
        the total number of identical iterations that will be run
    
    Returns
    -------
    n_workers_use : int
        the best choice for the number of workers
    n_worker_iters : int
        the maximum number of iterations allocated per worker at `n_workers_use`
    """

    # The maximum number of iterations per thread
    n_worker_iters = n_total // max_workers + (1 if n_total % max_workers > 0 else 0)
    # The maximum number of iterations per thread
    n_workers_use = int( np.ceil( n_total / n_worker_iters ) )
    
    return n_workers_use, n_worker_iters

def edges( xs: ArrayLike ) -> zip:
    """
    Iterate over paired left and right sides from lsit of bin edges `xs`
    """
    return zip( np.array( xs[:-1] ), np.array( xs[1:] ) )

def centers( xs: ArrayLike ) -> NDArray:
    """
    Centers of the bins with edges `xs`
    """
    return 0.5 * (np.array( xs[:-1] ) + np.array( xs[1:] ))

def bucket_list( window: 'tuple[float, float]',
                 divisions: int ) -> 'list[float]':
    """
    Bin edges for `window` divided into `divisions` different bins
    """
    ret = []
    spacing = (window[1] - window[0]) / divisions
    for i in range( divisions ):
        ret.append( (window[0] + i * spacing, window[0] + (i+1) * spacing) )
    return ret

def overlap(
        x: DataFrame,
        y: DataFrame,
        on: str = 'grid_dataset_20',
        key: str = 'ratio_active'
    ) -> 'tuple[int, int, int]':
    """
    Compute statistics about the overlap between the value of `key` in `x` and `y`, joined by `on`

    Parameters
    ----------
    x, y : DataFrame
        the input datasets
    on : str
        the key to use to align data points in `x` and `y` (default:
        'grid_dataset_20')
    key : str
        the key in `x` and `y` to look for overlapping nonzero values in
        (default: 'ratio_active')
    
    Returns
    -------
    overlap : int
        the number of values of `key` that are > 0 for both `x` and `y`
    either : int
        the number of values of `key` that are > 0 for either `x` or `y`
    possible : int
        the total number of values for `on` present in either `x` or `y` (i.e.,
        through an outer join on `on`)
    """
    
    keys = [on, key]
    
    combined = pd.merge( x[keys], y[keys],
                         on = on,
                         how = 'outer',
                         suffixes = ('_x', '_y') )
    
    possible = combined.shape[0]
    x_bool = ~np.isnan( combined[f'{key}_x'] ) & (combined[f'{key}_x'] > 0.)
    y_bool = ~np.isnan( combined[f'{key}_y'] ) & (combined[f'{key}_y'] > 0.)
    either = np.sum( x_bool | y_bool )
    overlap = np.sum( x_bool & y_bool )
    
    return overlap, either, possible

def _render_yaml_path(
        ps: Union[Path, str, 'list[str]']
    ) -> Path:
    """
    Normalize `ps`, which could be a `Path`, string, or list of path components, into a `Path`

    Raises
    ------
    ValueError
        if `ps` is an empty list
    """

    if type( ps ) == Path:
        return ps

    if type( ps ) == str:
        return Path( ps )

    # `ps` must be a list[str] here

    if len( ps ) < 1:
        raise ValueError( 'No path components specified' )

    ret = Path( ps[0] )
    for p in ps[1:]:
        ret = ret / Path( p )
    return ret

def _normalize_params( params: dict ) -> dict:
    """
    Perform default normalization of specific parameters loaded into `params`

    Returns a copy.
    """

    ret = params.copy()

    # Standard directories
    if 'helper_configs' in params:
        ret['helper_configs'] = [ _render_yaml_path( p )
                                  for p in params['helper_configs'] ]
    if 'hive_root' in params:
        ret['hive_root'] = _render_yaml_path( params['hive_root'] )
    if 'output_parent' in params:
        ret['output_parent'] = _render_yaml_path( params['output_parent'] )

    # Standard windows
    if 'responder_params' in params:
        if 'window_pre' in params['responder_params']:
            ret['responder_params']['window_pre'] = tuple( params['responder_params']['window_pre'] )
        if 'window_post' in params['responder_params']:
            ret['responder_params']['window_post'] = tuple( params['responder_params']['window_post'] )

    return ret

def load_params( p: Union[str, Path] ) -> dict:
    """
    Load the parameters from a `yaml` config file at path `p`, with extra spices

    Parameters
    ----------
    p : Union[str, Path]
        path to a `yaml` file, rendered as a string or `Path` object
    
    Returns
    -------
    dict
        the contents of the loaded `yaml` document at `p`, with some extra
        normalization steps handled for convenience
    """

    # Normalization
    if type( p ) == str:
        p = Path( p )
    
    with open( p, 'r' ) as f:
        ret = yaml.safe_load( f )
    ret = _normalize_params( ret )

    return ret