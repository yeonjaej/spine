"""Class to build all require representations."""

from collections import defaultdict, OrderedDict

import numpy as np

from spine.utils.logger import logger
from spine.utils.globals import COORD_COLS, VALUE_COL

from .fragment import FragmentBuilder
from .particle import ParticleBuilder
from .interaction import InteractionBuilder


class BuildManager:
    """Manager which constructs data representations based on the chain output.

    Takes care of two scenarios:
      - Interpret the raw output of the reconstruction chain
      - Load up existing objects stored as dictionaries
    """
    # Name of input data products needed to build representations. These names
    # are not set in stone, so they can be set in the configuration
    sources = {
            'data_tensor': ['data_adapt', 'data'],
            'label_tensor': 'clust_label',
            'label_adapt_tensor': ['clust_label_adapt', 'clust_label'],
            'label_g4_tensor': 'clust_label_g4',
            'depositions_q_label': 'charge_label',
            'sources': ['sources_adapt', 'sources'],
            'sources_label': 'sources_label',
            'particles': 'particles',
            'neutrinos': 'neutrinos'
    }

    def __init__(self, fragments, particles, interactions,
                 mode='both', units='cm', sources=None):
        """Initializes the build manager.

        Parameters
        ----------
        fragments : bool
            Build/load RecoFragment/TruthFragment objects
        particles : bool
            Build/load RecoParticle/TruthParticle objects
        interactions : bool
            Build/load RecoInteraction/TruthInteraction objects
        mode : str, default 'both'
            Whether to construct reconstructed objects, true objects or both
        sources : Dict[str, str], optional
            Dictionary which maps the necessary data products onto a name
            in the input/output dictionary of the reconstruction chain.
        """
        # Check on the mode, store desired units
        assert mode in ['reco', 'truth', 'both', 'all'], (
                f"Run mode not recognized: {mode}. Must be one of 'reco', "
                 "'truth', 'both' or 'all'.")
        self.mode = mode
        self.units = units
        
        # Parse the build sources based on defaults
        if sources is not None:
            for key, value in sources.items():
                assert key in self.sources, (
                         "Unexpected data product specified in `sources`: "
                        f"{key}. Should be one of {list(self.sources.keys())}.")
            self.sources.update(**sources)

        for key, value in self.sources.items():
            if isinstance(value, str):
                self.sources[key] = [value]
            else:
                self.sources[key] = value

        # Initialize the builders
        self.builders = OrderedDict()
        if fragments:
            self.builders['fragment'] = FragmentBuilder(mode, units)
        if particles:
            self.builders['particle'] = ParticleBuilder(mode, units)
        if interactions:
            assert particles, (
                    "Interactions are built from particles. If "
                    "`interactions` is True, so must "
                    "`particles` be.")
            self.builders['interaction'] = InteractionBuilder(mode, units)

        assert len(self.builders), (
                "Do not call the builder unless it does anything.")

    def __call__(self, data):
        """Build the representations for one entry.

        Parameters
        ----------
        data : dict
            Dictionary of input data and model outputs

        Notes
        -----
        Modifies the data dictionary in place.
        """
        # If this is the first time the builders are called, build
        # the objects shared between fragments/particles/interactions
        load = True
        if 'points' not in data:
            load = False
            if np.isscalar(data['index']):
                sources = self.build_sources(data)
            else:
                sources = defaultdict(list)
                for entry in range(len(data['index'])):
                    sources_e = self.build_sources(data, entry)
                    for key, val in sources_e.items():
                        sources[key].append(val)

            data.update(**sources)

        # Loop over builders
        for name, builder in self.builders.items():
            # Build representations
            builder(data)

            # Generate match pairs from stored matches
            if load and self.mode in ['both', 'all']:
                if np.isscalar(data['index']):
                    match_dict = self.load_match_pairs(data, name)
                else:
                    match_dict = defaultdict(list)
                    for entry in range(len(data['index'])):
                        match_dict_e = self.load_match_pairs(data, name, entry)
                        for key, val in match_dict_e.items():
                            match_dict[key].append(val)

                data.update(**match_dict)

    def build_sources(self, data, entry=None):
        """Construct the reference coordinate and value tensors used by
        all the representations built by the module.

        These objects should be stored along with the constructed objects
        if the objects are to be loaded later on.

        Parameters
        ----------
        data : dict
            Dictionary of input data and model outputs
        entry : int, optional
            Entry number
        """
        # Fetch the orginal sources
        sources = {}
        for key, alt_keys in self.sources.items():
            for alt in alt_keys:
                if alt in data:
                    sources[key] = data[alt]
                    if entry is not None:
                        sources[key] = data[alt][entry]
                    break

        # Build aditional information
        update = {}
        update['points'] = sources['data_tensor'][:, COORD_COLS]
        update['depositions'] = sources['data_tensor'][:, VALUE_COL]
        
        if self.mode != 'reco':
            update['label_tensor'] = sources['label_tensor']
            update['points_label'] = sources['label_tensor'][:, COORD_COLS]
            update['depositions_label'] = sources['label_tensor'][:, VALUE_COL]

            update['label_adapt_tensor'] = sources['label_adapt_tensor']
            update['depositions_label_adapt'] = (
                    sources['label_adapt_tensor'][:, VALUE_COL])

            if 'depositions_q_label' in sources:
                update['depositions_q_label'] = (
                        sources['depositions_q_label'][:, VALUE_COL])

        if 'label_g4_tensor' in sources:
            update['label_g4_tensor'] = sources['label_g4_tensor']
            update['points_g4'] = sources['label_g4_tensor'][:, COORD_COLS]
            update['depositions_g4'] = sources['label_g4_tensor'][:, VALUE_COL]

        if 'sources' in sources:
            update['sources'] = sources['sources'].astype(int)
        if 'sources_label' in sources:
            update['sources_label'] = sources['sources_label'].astype(int)

        # If provided, etch the point attributes to check their units
        for obj in ['fragment', 'particle']:
            for key in [f'{obj}_start_points', f'{obj}_end_points']:
                if key in data:
                    update[key] = data[key]
                    if entry is not None:
                        update[key] = update[key][entry]

        # Convert everything to the proper units once and for all
        if self.units != 'px':
            # Fetch metadata
            assert 'meta' in data, (
                    "Must provide metadata to build objects in cm.")

            meta = data['meta'][entry] if entry is not None else data['meta']
            for key in update:
                if 'points' in key and key in update:
                    if key in update:
                        update[key] = meta.to_cm(
                                np.copy(update[key]), center=True)

            for key in ['particles', 'neutrinos']:
                if key in sources:
                    update[key] = sources[key]
                    for obj in sources[key]:
                         if obj.units != self.units:
                             obj.to_cm(meta)

        return update

    @staticmethod
    def load_match_pairs(data, name, entry=None):
        """Generate lists of matched object pairs from stored matches.

        Parameters
        ----------
        data : dict
            Dictionary of input data and model outputs
        name : str
            Object type name
        entry : int, optional
            Entry number
        """
        # Initialize the name of the match lists
        prefix = f'{name}_matches'

        # Create match pairs in both directions (true to reco and vice versa)
        result = {}
        for source, target in [('reco', 'truth'), ('truth', 'reco')]:
            # Fetch the lists of objects to match
            sources = data[f'{source}_{name}s']
            targets = data[f'{target}_{name}s']
            if entry is not None:
                sources, targets = sources[entry], targets[entry]

            # Loop
            suffix = f'{source[0]}2{target[0]}'
            match_key = f'{prefix}_{suffix}'
            match_overlap_key = f'{match_key}_overlap'
            result[match_key] = []
            result[match_overlap_key] = []
            for obj in sources:
                if not obj.is_matched:
                    # If no match is found, give an empty value to the match
                    result[match_key].append((obj, None))
                    result[match_overlap_key].append(-1.)

                else:
                    # If a match is found, the first is always the best match
                    best_match = obj.match_ids[0]
                    result[match_key].append((obj, targets[best_match]))
                    result[match_overlap_key].append(obj.match_overlaps[0])

        return result
