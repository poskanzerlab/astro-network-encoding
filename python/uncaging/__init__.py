"""
Tools for interfacing with AQuA datasets from uncaging experiments

Designed to accompany the paper "Network-level encoding of local
neurotransmitters in cortical astrocytes" (DOI: TODO)
"""


## ------------------------------------------------------------------------ ##
## Imports
## ------------------------------------------------------------------------ ##

# Type hints
from typing import (
    Any,
    Tuple,
    Optional,
    Union,
    Callable,
    Protocol,
)

from numpy.typing import (
    ArrayLike,
)

from matplotlib.pyplot import (
    Axes,
)

from aqua.manager import (
    HiveManager,
)

# Standard library
import os
import time
import itertools

# Packages
import yaml
from tqdm import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.stats
import scipy.linalg as sla
import scipy.integrate as sint

import aqua.manager as aman
import aqua.stats as astats

# Other uncaging components
from . import utils
from . import events

# One-offs
from functools import reduce
from pathlib import Path
from datetime import datetime

from pandas import DataFrame
from numpy import ndarray


## ------------------------------------------------------------------------ ##
## Protocols (for typing)
## ------------------------------------------------------------------------ ##

class DFStatistic( Protocol ):
    "A statistic function, which produces a `float` from a `DataFrame`"
    def __call__(
            self,
            df: DataFrame
        ) -> float:
        return super().__call__( df )

class FigureSaver( Protocol ):
    "A function that saves a figure with the desired `title` and `extension`"
    def __call__(
            self,
            title: str,
            extension: str
        ) -> None:
        return super().__call__( title, extension )

class DataSaver( Protocol ):
    "A function that saves `numpy`-formatted data with the given `title` and args"
    def __call__(
            self,
            title: str,
            *args,
            **kwargs
        ) -> None:
        return super().__call__( title, *args, **kwargs )

class CSVSaver( Protocol ):
    "A function that saves out a `DataFrame` to a `csv` with the specified `title`"
    def __call__(
            self,
            title: str,
            data: DataFrame
        ) -> None:
        return super().__call__( title, data )


## ------------------------------------------------------------------------ ##
## Functions
## ------------------------------------------------------------------------ ##

## ** Private ** ##

## Filename spec formatting

def _unwrap_yaml_filename_spec( raw_spec: list ) -> list:
    """
    Converts the YAML filename specification format to the Python equivalent

    YAML only supports lists (and not tuples); this takes each terminal entry of
    the spec (as lists of lists) and appropriately converts them to tuples

    Parameters
    ----------
    raw_spec : list
        a filename specification in the YAML-supported nested list format

    Returns
    -------
    list
        a filename specification in the native `uncaging.py` format
    
    Raises
    ------
    ValueError
        if the `raw_spec` is somehow malformed
    """

    ret = []

    for item in raw_spec:
        if type( item ) == list:
            if type( item[0] ) == str:
                # Terminal entry: `item` is specifying a single part
                new_item = tuple( item )
            elif type( item[0] ) == list:
                # Entry contains subitems: must recurse
                new_item = _unwrap_yaml_filename_spec( item )

            ret.append( new_item )
        
        else:
            raise ValueError( 'Each raw YAML spec entry must be a list' )

    return ret

ReplacementSpec = 'dict[str, Any]'

def _replace_objects(
        t: tuple,
        **replacements: ReplacementSpec,
    ) -> tuple:
    """
    A second method for doing parameter replacement in format spec tuples

    TODO: There is a better unified way to do dynamic options in `yaml` config
    files, this hack was a quick compromise

    Parameters
    ----------
    t : tuple
        a tuple of inputs; the strings in it will be altered
    **replacements : ReplacementSpec
        each instance of '$k' in an entry of `t` (for 'k' a key in
        `replacements`) will be replaced with `replacements[k], the
        corresponding value
    
    Returns
    -------
    tuple
        the same format as `t`, but with the substitutions made
    """
    return tuple( ( replacements[ti[1:]]
                    if type( ti ) == str and ti[0] == '$'
                    else ti )
                  for ti in t )

def _format_strings(
        t: tuple,
        **replacements: ReplacementSpec,
    ) -> tuple:
    """
    Take any format string in `t` and perform the string substitution in `replacements`
    
    Parameters
    ----------
    t : tuple
        a tuple of inputs; the strings in it will be formatted
    **replacements : dict
        each instance of '{k}' in an entry of `t` (for 'k' a key in
        `replacements`) will be replaced with `replacements[k]`, the
        corresponding value

    Returns
    -------
    tuple
        the same format as `t`, but with the substitutions made
    """
    return tuple( ( ti.format( **replacements )
                    if type( ti ) == str
                    else ti )
                  for ti in t )

def _spec_replace(
        spec: list,
        **replacements: ReplacementSpec,
    ) -> list:
    """
    Recursively reformat the entries specification `spec`
    
    Parameters
    ----------
    spec : list
        a full `filename_spec` object: a list whose entries are either
        * a tuple (giving one specification for one filename component), or
        * a list (of multiple specifications that are all run on the
            same filename component)
    **replacements : dict
        each instance of '{k}' in an entry or subentry of `spec` (for 'k' a key
        in `replacements`) will be replaced with `replacements[k]`, the
        corresponding value

    Returns
    -------
    list
        the same format as `spec`, but with the substitutions made
    """
    ret = [ ( _replace_objects( _format_strings( x, **replacements ), **replacements )
              if type( x ) == tuple
              else _spec_replace( x, **replacements ) )
            for x in spec ]
    return ret


## ------------------------------------------------------------------------ ##
## `Helper`` class
## ------------------------------------------------------------------------ ##

# Some type shortcuts, for brevity
HelperConfigPath = 'Optional[Union[str, Path, list[Union[str, Path]]]]'
HelperOutputPath = Optional[Union[str, Path]]
HelperFilenameSpec = 'list[Union[tuple, list[tuple]]]'
HelperManagerSpec = Union['dict[str, HiveManager]', 'list[HiveManager]']
HelperLoadEventsResult = Union['dict[str, DataFrame]', 'list[DataFrame]']

##

class Helper:
    """
    Encapsulates information about one experiment's ins, outs, and configuration
    """

    def __init__(
            self,
            config: Optional[dict] = None,
            config_path: HelperConfigPath = None,
            hive_root: str = '.',
            output_to: HelperOutputPath = None,
        ) -> None:
        """
        Load the configuration for a new analysis helper.

        Parameters
        ----------
        config : dict, optional
            the configuration data to load for the analysis; if `None`
            (default), looks to `config_path` to find where to load the config
            data from
        config_path : HelperConfigPath
            either:
            * a path to a dir or `.yaml` file to load config from; or,
            * a list of these.
            If not specified, must manually provide a `config` object.
        hive_root : str, optional
            path to the root of the input data hive (default: '.', the current
            working directory)
        output_to : str, optional
            path to the parent directory to place outputs in; if not provided,
            output-saving functions will raise exceptions

        Raises
        ------
        ValueError
            if no `config` or `config_path` is specified, or if `config_path` is
            poorly specified
        """

        if config is None:
            if config_path is None:
                raise ValueError( 'Must specify a `config_path` or a `config` object' )
            
            if type( config_path ) != list:
                config_path = [config_path]

            # `config_path` normalized as a list
            config = dict()
            for cur_path in config_path:

                if type( cur_path ) == str:
                    cur_path = Path( cur_path )
                
                if not issubclass( type( cur_path ), Path ):
                    raise ValueError( f'Unsupported type for a config path: {type( cur_path )}' )

                # `cur_path` normalized as a Path
                cur_contents = []
                if cur_path.is_dir():
                    cur_config = dict()
                    # Iterate through dir contents
                    for p in cur_path.glob( '*.yaml' ):
                        cur_contents.append( p )
                else:
                    cur_contents.append( cur_path )
                
                # `cur_contents` normalized as list of paths
                cur_config = dict()
                for p in cur_contents:
                    with open( p, 'r' ) as f:
                        cur_config = utils._merge_dicts( cur_config, yaml.safe_load( f ) )
                
                # Now we have `cur_config` fully; integrate it
                config = utils._merge_dicts( config, cur_config )
        
        self._config = config

        if output_to is None:
            self._output_parent = None
        else:
            if issubclass( type( output_to ), Path ):
                self._output_parent = output_to
            elif type( output_to ) == str:
                self._output_parent = Path( output_to )
            else:
                raise ValueError( f'Unsupported type for output directory: {type( output_to )}' )
        
        self._hive_root = Path( hive_root )

    @property
    def condition_labels( self ) -> 'Optional[dict[str, str]]':
        """
        TODO Is this used externally?
        """
        return self._config.get( 'condition_labels', None )

    @property
    def condition_colors( self ) -> 'dict[str, str]':
        """
        `matplotlib` color specification for each experimental condition
        """

        if 'condition_colors' not in self._config:
            # TODO Refine Exception
            raise Exception( "'conditon_colors' not specified in Helper config" )

        return self._config['condition_colors']

    @property
    def mark_descriptions( self ) -> 'dict[str, str]':
        """
        Human-readable description for each mark
        """

        if 'mark_descriptions' not in self._config:
            # TODO Refine Exception
            raise Exception( "'mark_descriptions' not specified in Helper config" )

        return self._config['mark_descriptions']

    def _experiment_spec(
            self,
            experiment: str
        ) -> dict:
        """
        Extract the config for a particular experiment from the imported configs

        Parameters
        ----------
        experiment : str
            The name of the experiment to query
        
        Returns
        -------
        dict
            Config data for the experiment, loaded from all of this `Helper`'s
            `yaml` files
        
        Raises
        ------
        KeyError
            no config is provided for the specified `experiment`
        Exception
            multiple matching configs provided for the specified `experiment`
            TODO refine this exception
        """

        matching_experiments = [ e
                                 for e in self._config['experiments']
                                 if e['name'] == experiment ]
        
        if len( matching_experiments ) < 1:
            raise KeyError( f"No Helper config provided for experiment named '{experiment}'" )

        if len( matching_experiments ) > 1:
            # TODO Refine
            raise Exception( f"Multiple matching configs provided for experiment '{experiment}'" )
        
        return matching_experiments[0]

    @property
    def _filename_spec_template( self ) -> HelperFilenameSpec:
        """
        Generic template components for parsing filenames in a hive

        Shortcut for calling `Helper._filename_spec` with no arguments

        Returns
        -------
        spec_template : HelperFilenameSpec
            see `Helper._filename_spec`
        """
        spec_template = self._filename_spec()
        return spec_template

    def _filename_spec(
            self,
            experiment: Optional[str] = None
        ) -> HelperFilenameSpec:
        """
        List of components for parsing filenames in a hive

        Parameters
        ----------
        experiment : str, optional
            if not specified, will return the generic spec template; if
            specified, will use the experiment configuration to do substitutions
            on this template specific to that experiment's data
        
        Returns
        -------
        HelperFilenameSpec
            the filename specification. Each entry is either
            * a tuple describing how this part of the filename should be
            parsed; or,
            * a list of these tuples, each of which is executed on the same
            part
        
        """

        spec_unwrapped = _unwrap_yaml_filename_spec( self._config['filename_spec'] )
        if experiment is None:
            return spec_unwrapped
        
        # Experiment is specified, so can perform specific replacement
        experiment_spec = self._experiment_spec( experiment )
        spec_replaced = _spec_replace( spec_unwrapped,
                                       condition_labels = self.condition_labels,
                                       suffix = experiment_spec['hive_suffix'] )
        return spec_replaced
    
    def _output_dir(
            self,
            *descriptors: Tuple[str, ...],
            subpath: Optional[Path] = None,
            timestamp: Optional[str] = None
        ) -> Path:
        """
        Get the path for an analysis output directory

        Parameters
        ----------
        *descriptors : Tuple[str, ...]
            additional descriptor strings to append to the output directory name
        subdir : Path, optional
            the subdirectory to create under the main output parent directory
        timestmap : str, optional
            a string timestmap representation to avoid collisions; default is to
            use the current date/time as %Y-%m-%d-%H-%M-%S

        Returns
        -------
        Path
            the path to the requested output dir
        
        Raises
        ------
        FileNotFoundError
            if no output parent dir is specified in tihs `Helper`
        """

        if self._output_parent is None:
            raise FileNotFoundError( 'No output parent dir specified' )
        
        if timestamp is None:
            now = datetime.now()
            timestamp = now.strftime( '%Y_%m_%d_%H_%M_%S' )
        
        components = [a for a in descriptors] + [timestamp]
        dirname = '-'.join( components )
        
        if subpath is None:
            return self._output_parent / dirname
        return self._output_parent / subpath / dirname
    
    def _make_output_dir(
            self,
            *descriptors: Tuple[str, ...],
            subpath: Optional[Path] = None,
            timestamp: str = None
        ) -> Path:
        """
        Make (and return) a directory for storing analysis outputs

        Parameters
        ----------
        *descriptors : tuple[str, ...]
            additional descriptor strings to append to the output directory name
        subdir : Path, optional
            the subdirectory to create under the main output parent directory
        timestmap : str, optional
            a string timestmap representation to avoid collisions; default is to
            use the current date/time as %Y-%m-%d-%H-%M-%S

        Returns
        -------
        Path
            the path to the requested output dir
        """
        
        output_dir = self._output_dir( *descriptors,
                                       subpath = subpath,
                                       timestamp = timestamp )
        os.makedirs( output_dir,
                     exist_ok = True )
        return output_dir
    
    def figure_saver(
            self,
            stem: str,
            analysis: str
        ) -> FigureSaver:
        """
        Make a figure-saving function for the current analysis

        Parameters
        ----------
        stem : str
            filename stem to add to the beginning of every saved output
        analysis: str
            name of the particular analysis to add to the end of every saved output

        Returns
        -------
        FigureSaver
            the figure-saving function
        """

        output_dir = self._make_output_dir( stem, analysis,
                                            subpath = Path( 'figures' ) )
        
        def ret( title: str, extension: str = 'svg' ):
            """
            Save a figure named `title` with the given `extension`
            """
            components = [stem, title, analysis]
            filename = '_'.join( components ) + f'.{extension}'
            output_path = output_dir / filename
            plt.savefig( output_path )
        
        return ret
    
    def data_saver(
            self,
            stem: str,
            analysis: str
        ) -> DataSaver:
        """
        Make a data-saving function for the current analysis

        Parameters
        ----------
        stem : str
            filename stem to add to the beginning of every saved output
        analysis: str
            name of the particular analysis to add to the end of every saved output

        Returns
        -------
        DataSaver
            the data-saving function
        """

        output_dir = self._make_output_dir( stem, analysis,
                                            subpath = Path( 'arrays' ) )
        
        def ret( title: str, *args, **kwargs ):
            """
            Save data in `*args` or `**kwargs` to an `npz` file named `title`
            
            Parameters
            ----------
            title : str
                the title for the newly saved data file
            *args : tuple
                passed to `np.savez`
            **kwargs : dict
                passed to `np.savez`
            """
            components = [stem, title, analysis]
            filename = '_'.join( components ) + '.npz'
            output_path = output_dir / filename
            np.savez( output_path, *args, **kwargs )
        
        return ret

    def csv_saver(
            self,
            stem: str,
            analysis: str,
        ) -> CSVSaver:
        """
        Make a `pandas` csv-saving function for the current analysis

        Parameters
        ----------
        output_dir : str
            output figure directory
        stem : str
            filename stem to add to the beginning of every saved output
        analysis: str
            name of the particular analysis to add to the end of every saved output

        Returns
        -------
        CSVSaver
            the csv-saving function
        """

        output_dir = self._make_output_dir( stem, analysis,
                                            subpath = Path( 'tables' ) )
        
        def ret( title: str, data: DataFrame ):
            """
            Save `data` to a `.csv` file named `title`
            
            Parameters
            ----------
            title : str
                the title for the newly saved data file
            data : pandas.DataFrame
                data to be saved out
            """
            components = [stem, title, analysis]
            filename = '_'.join( components ) + '.csv'
            output_path = os.path.join( output_dir, filename )
            data.to_csv( output_path )
        
        return ret

    def _get_manager(
            self,
            experiment: str,
            labels: Optional[Any] = None,
        ) -> HiveManager:
        """
        Use the loaded config to make a `HiveManager` for the data for the given `experiment`

        Parameters
        ----------
        experiment : str
            The name of the experiment to generate a manager for
        labels
            TODO

        Returns
        -------
        aqua.manager.HiveManager
            An `aqua` object for managing the events for `experiment`
        """

        ## Get the manager

        experiment_spec = self._experiment_spec( experiment )

        # Set up the filename parser for managing our data
        filename_spec = self._filename_spec( experiment = experiment )
        parser = aman.get_standard_filename_parser( filename_spec )
        
        # Join together the path components for the experiment's relative path
        experiment_rel_path = Path( os.path.join( *tuple( experiment_spec['hive_path'] ) ) )

        manager = aman.HiveManager( self._hive_root / experiment_rel_path,
                                    parser = parser )

        ## Add additional spices

        # Label each dataset with the proper orientation as specified in a file
        if labels is not None:
            directionality_labels = aman.load_labels( labels,
                                                      key = 'left_right' )
            manager.add_dataset_postprocessor( 'left_right',
                                               aman.get_label_postprocessor( directionality_labels ) )
            manager.add_laterality_postprocessors()

        return manager

    def _get_managers(
            self,
            experiments: 'Optional[list[str]]' = None,
            comparisons: 'Optional[list[tuple[str, str]]]' = None
        ) -> 'dict[str, HiveManager]':
        """
        Make and collate `HiveManager`s for the specified `experiments` or `comparisons`

        Parameters
        ----------
        experiments : list[str], optional
            a list of experiments laid out in the Helper config
        comparisons : list[tuple[str, str]], optional
            a list of tuples of the form `(a, b)` where each of `a` and `b` is
            an experiment (for when the calling script is using a list of
            comparisons, for example)
        
        Returns
        -------
        dict[str, aqua.manager.HiveManager]
            the `HiveManager` for each experiment implied by either
            `experiments` or `comparisons`
        """
        
        # Infer the experiments that would be implied by `comparisons`
        if comparisons is not None:
            # ... which are the unique experiments given in `comparisons`
            comparison_experiments_unique = set( reduce( lambda a, b: a + b,
                                                         comparisons,
                                                         tuple() ) )
        else:
            comparison_experiments_unique = set()

        # Infer all of the available experiments
        if experiments is None:
            experiments_unique = set( [ e['name']
                                        for e in self._config['experiments'] ] )
        else:
            experiments_unique = set( experiments )
        
        # Take all of the experiments implied by either condition
        experiments_all = list( experiments_unique
                                | comparison_experiments_unique )

        ret = { e: self._get_manager( e )
                for e in experiments_all }
        
        return ret

    def _load_events(
            self,
            managers: HelperManagerSpec,
            #   use_experiments: Union[list, bool] = True,
            **kwargs
        ) -> HelperLoadEventsResult:
        """
        TODO

        Parameters
        ----------
        managers : HelperManagerSpec
            TODO
        
        ...
        
        Returns
        -------
        HelperLoadEventsResult
            TODO
        """
        
        if type( managers ) == dict:

            ret_headers = dict()
            ret_events = dict()

            for k, manager in managers.items():
                
                ret_headers[k], ret_events[k] = manager.all_events( **kwargs )
                
                cur_experiment_hits = [ e
                                        for e in self._config['experiments']
                                        if e['name'] == k ]
                cur_experiment = cur_experiment_hits[0]
                
                query = np.array( [ True
                                    for _ in range( ret_events[k].shape[0] ) ] )
                
                for qk, qv in cur_experiment['query']:
                    cur_query = (ret_events[k][qk] == qv)
                    query = query & cur_query
                
                ret_events[k] = ret_events[k][query]
        
        elif type( managers ) == list:
            raise NotImplemented()
                    
        # TODO Is this used? It doesn't seem referenced
        #     ret_headers = []
        #     ret_events = []
            
        #     for i_manager, manager in enumerate( managers ):
        #         h, e = manager.all_events( **kwargs )
                
        #         if type( use_experiments ) == list:
        #             queries = experiments[i_manager]['query']
        #             filter_query = np.array( [True for _ in ret_events[k]] )
        #             for column, value in queries:
        #                 filter_query = filter_query & (ret_events[i_manager][ret_events[i_manager][column] == value] )
                    
        #             e = e[filter_query]
                
        #         ret_headers.append( h )
        #         ret_events.append( e )
        
        else:
            raise TypeError( f'Invalid type for managers: {type( managers )}' )
        
        return ret_headers, ret_events

    def load_analysis_events(
            self,
            analysis: str,
            verbose: bool = False,
            extra_decorators: 'list[tuple]' = None,
        ) -> DataFrame:
        """
        Use `aqua-py` to load events for a specific analysis

        Parameters
        ----------
        analysis : {'wt', 'cx43', 'invivo', 'wt-repeat', 'wt-repeat-70'}
            the analysis to load
        verbose : bool, optional
            if `True`, print additional debug messages (default: False)
        extra_decorators : list[tuple], optional
            a list of additional `aqua-py` event decorators to compute when
            extracting the data; see TODO (default is a 'grid' with a scale that is
            reasonable for the specified `analysis`; set to `[]` if you don't need
            the grid, which can be time-consuming to compute)

        Returns
        -------
        pandas.DataFrame
            all of the events for the specified `analysis`
        
        Raises
        ------
        ValueError
            if there are multiple experiments in the config that have the same
            `name`, but specify different templates for ramping cells
        """
        
        # Parse config
        # ------------

        config_post_defaults = {
            'dataset_keys': [],
            'header_keys': [],
            'decorators': [],
            'extra_decorators_default': [],
            'exclude_ramp': False,
            'ramp_group_key': 'cell_global',
            'ramp_n_events_thresh': 5,
            'ramp_p_thresh': 0.1,
        }

        config_post = self._config['postprocessing']
        for k, v in config_post_defaults.items():
            if k not in config_post:
                config_post[k] = v
        
        dataset_keys = config_post['dataset_keys']
        header_keys = config_post['header_keys']

        default_decorators = [ tuple( d )
                               for d in config_post['decorators'] ]
        if extra_decorators is None:
            extra_decorators = [ tuple( d )
                                for d in config_post['extra_decorators_default'] ]
        decorators = default_decorators + extra_decorators

        global_decorators = [ tuple( d )
                              for d in config_post['global_decorators'] ]
        
        # Normalization
        for d in decorators:
            if type( d ) != tuple:
                d = tuple( d )
        for d in global_decorators:
            if type( d ) != tuple:
                d = tuple( d )

        # Determine if we need to add coregistered grid decorators
        i_coreg_grids = None
        coreg_grids_by = None
        for i, gd in enumerate( global_decorators ):
            if gd[0] == 'coreg_grids':
                i_coreg_grids = i
                coreg_grids_by = gd[1]
                break
        
        if coreg_grids_by is not None:
            # Get rid of the coreg grids decorator for normalization
            del global_decorators[i_coreg_grids]

            # Add a matching coreg decorator for each grid
            for d in decorators:
                if d[0] == 'grid':
                    grid_key_local = f'grid_dataset_{d[1]:0d}'
                    grid_key_coreg = f'grid_global_{d[1]:0d}_coreg'
                    new_decorator = ('coreg', grid_key_coreg, coreg_grids_by + [grid_key_local])
                    global_decorators.append( new_decorator )

        experiment_names = self._config['analyses'][analysis]['experiments']
        experiment_managers = self._get_managers( experiments = experiment_names )
        # This ensures that `experiment_names` actually reflects what we have
        experiment_names = list( experiment_managers.keys() )

        # Load events for each experiment
        # -------------------------------
        experiment_headers, experiment_events = self._load_events(
            experiment_managers,
            verbose = verbose,
            dataset_keys = dataset_keys,
            header_keys = header_keys,
        )

        # Add single-dataset decorators
        # -----------------------------
        experiment_events = { k: events._add_event_decorators(
                                    experiment_headers[k],
                                    experiment_events[k],
                                    decorators = decorators
                                )
                              for k in experiment_events.keys() }

        # Determine ramping cells
        # -----------------------
        if config_post['exclude_ramp']:

            df_ramp_all = None

            for experiment_name in experiment_names:
                
                cur_experiment_hits = [ e
                                        for e in self._config['experiments']
                                        if e['name'] == experiment_name ]
                cur_ramp_path_templates = [ os.path.join( *tuple( e['ramp_path_template'] ) )
                                            for e in cur_experiment_hits ]
                
                ramp_path_template = cur_ramp_path_templates[0]
                for x in cur_ramp_path_templates[1:]:
                    if x != ramp_path_template:
                        raise ValueError( f'Non-matching ramp path templates for experiment: {experiment_name}' )
                
                ramp_path_cur = ramp_path_template.format( experiment = experiment_name )
                df_ramp = pd.read_csv( ramp_path_cur )
                df_ramp['experiment'] = experiment_name

                if df_ramp_all is None:
                    df_ramp_all = df_ramp
                else:
                    df_ramp_all = pd.concat( [df_ramp_all, df_ramp] )

            # Exclude events from ramping cells in each experiment, if needed
            # --------------------------------------------
            events_all = None

            for experiment_name in experiment_names:

                groups_exclude = []
                groups_include = []

                filter_ramp_spec_experiment = df_ramp_all['experiment'] == experiment_name
                df_ramp_cur = df_ramp_all[filter_ramp_spec_experiment]

                # Separate out the event groupings that should be excluded or
                # included on the basis of the ramping criteria
                for _, row in df_ramp_cur.iterrows():
                    if ( row['n_events'] <= config_post['ramp_n_events_thresh']
                         or row['p'] < config_post['ramp_p_thresh'] ):
                        groups_exclude.append( row[config_post['ramp_group_key']] )
                    else:
                        groups_include.append( row[config_post['ramp_group_key']] )

                events_cur = experiment_events[experiment_name]
                # TODO This should be equivalent to using `groups_include`
                # rather than the inverse of `groups_exclude` membership, but
                # I'm going to leave it in the original form for consistency
                # with the manuscript execution
                filter_events_not_excluded = ~events_cur[config_post['ramp_group_key']].isin( groups_exclude )
                events_cur_filtered = events_cur[filter_events_not_excluded]

                if events_all is None:
                    events_all = events_cur_filtered
                else:
                    events_all = pd.concat( [events_all, events_cur_filtered],
                                            ignore_index = True )

                if verbose:
                    print()
                    print( f'Experiment {experiment_name}' )
                    print( f'    Excluded cells: {len( groups_exclude )}' )
                    print( f'    Included cells: {len( groups_include )}' )
                    print()
            
            #
            
        else:
            
            # Don't need to do any exclusions, so just concatenate everything
            events_all = pd.concat( [ experiment_events[e]
                                      for e in experiment_names ],
                                    ignore_index = True )
        
        # Add global decorators
        # ---------------------
        events_all = events._add_event_decorators(
            None,
            events_all,
            decorators = global_decorators
        )
        
        return events_all.copy()

#