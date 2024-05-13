"""Manages the operation of post-processors."""

import numpy as np
from warnings import warn
from collections import defaultdict, OrderedDict

from .factories import post_processor_factory

from mlreco.utils.stopwatch import StopwatchManager


class PostManager:
    """Manager in charge of handling post-processing scripts.
    
    It loads all the post-processor objects once and feeds them data.
    """

    def __init__(self, cfg, parent_path=''):
        """Initialize the post-processing manager.

        Parameters
        ----------
        cfg : dict
            Post-processor configurations
        parent_path : str, optional
            Path to the analysis tools configuration file
        """
        # Loop over the post-processor modules and get their priorities
        keys = np.array(list(cfg.keys()))
        priorities = -np.ones(len(keys), dtype=np.int32)
        for i, k in enumerate(keys):
            if 'priority' in cfg[k]:
                priorities[i] = cfg[k].pop('priority')

        # Add the modules to a processor list in decreasing order of priority
        self.watch = StopwatchManager()
        self.modules = OrderedDict()
        keys = keys[np.argsort(-priorities)]
        for k in keys:
            # Profile the module
            self.watch.initialize(k)

            # Append
            self.modules[k] = post_processor_factory(k, cfg[k], parent_path)

    def __call__(self, data):
        """Pass one batch of data through the post-processors.

        Parameters
        ----------
        data : dict
            Dictionary of data products
        """
        # Loop over the post-processor modules
        for key, module in self.modules.items():
            # Run the post-processor on each entry
            self.watch.start(key)
            if np.isscalar(data['index']):
                result = module(data)
            else:
                num_entries = len(data['index'])
                result = defaultdict(list)
                for entry in range(num_entries):
                    result_e = module(data, entry)
                    for key, value in result_e.items():
                       result[key].append(value)
            self.watch.stop(key)

            # Update the input dictionary
            for key, val in result.items():
                assert len(val) == num_entries
                data[key] = val
