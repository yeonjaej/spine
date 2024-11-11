"""Analysis module template.

Use this template as a basis to build your own analysis script. An analysis
script takes the output of the reconstruction and the post-processors and
performs basic selection cuts and store the output to a CSV file.
"""

# Add the imports specific to this module here
import numpy as np
import pandas as pd
import yaml, os, sys, re

# Must import the analysis script base class
from spine.ana.base import AnaBase

# Must list the post-processor(s) here to be found by the factory.
# You must also add it to the list of imported modules in the
# `spine.ana.factories`!
__all__ = ['pi0Ana']


class pi0Ana(AnaBase):
    """Script to make numu CC pi0 selection"""
    name = 'pi0'

    def __init__(self, **kwargs):
        """Initialize the analysis script.
        """
        # Initialize the parent class
        super().__init__('interaction', 'both', **kwargs)

        # Initialize the CSV writer(s) you want
        self.initialize_writer('log')
        
        self.keys['interaction_matches_r2t'] = True

    def process(self, data):
        """Pass data products corresponding to one entry through the analysis.

        Parameters
        ----------
        data : dict
            Dictionary of data products
        """
        
        # Loop over matched interactions (r2t)      
        interaction_matches_r2t = data['interaction_matches_r2t']
        for match in interaction_matches_r2t:
            
            # Get match components
            reco_inter = match[0]
            true_inter = match[1]
            if true_inter == None:
                continue
            
            # Containment cut
            if not reco_inter.is_contained : continue
            
            # Fiducial cut
            if not reco_inter.is_fiducial : continue
            
            # Flash cut
            if reco_inter.flash_time < 0 or reco_inter.flash_time > 1.6 : continue
            
            # Primary muon cut
            reco_muons = [p for p in reco_inter.particles if (p.pid == 2) and (p.is_primary)]
            if len(reco_muons) != 1 : continue
            
            # Primary photons cut
            reco_photons = [p for p in reco_inter.particles if (p.pid == 0) and (p.is_primary)]
            if len(reco_photons) < 2 : continue
            
            # Sort photons by energy
            reco_photons = sorted([rp for rp in reco_photons], key=lambda rp : rp.calo_ke, reverse=True)
            
            # This is our pi0 event
            reco_pi0_photons = reco_photons[:2]
            
            # Opening angle calculation and invariant mass
            # (use interaction vertex...only for CC events)
            '''
            reco_pi0_leading_photon_dir_s2v = (reco_pi0_photons[0].start_point - reco_inter.vertex) / np.linalg.norm(reco_pi0_photons[0].start_point - reco_inter.vertex)
            reco_pi0_subleading_photon_dir_s2v = (reco_pi0_photons[1].start_point - reco_inter.vertex) / np.linalg.norm(reco_pi0_photons[1].start_point - reco_inter.vertex)
            reco_pi0_cos_opening_angle_s2v = np.dot(reco_pi0_leading_photon_dir_s2v, reco_pi0_subleading_photon_dir_s2v)
            reco_pi0_opening_angle_s2v = np.degrees(np.arccos(reco_pi0_cos_opening_angle_s2v))
            reco_pi0_mass_s2v = np.sqrt(2*reco_pi0_photons[0].calo_ke*reco_pi0_photons[1].calo_ke*(1-reco_pi0_cos_opening_angle_s2v))
            '''
            
            # (use shower dirs...for CC or NC events)
            reco_pi0_leading_photon_dir = reco_pi0_photons[0].start_dir / np.linalg.norm(reco_pi0_photons[0].start_dir)
            reco_pi0_subleading_photon_dir = reco_pi0_photons[1].start_dir / np.linalg.norm(reco_pi0_photons[1].start_dir)
            reco_pi0_cos_opening_angle = np.dot(reco_pi0_leading_photon_dir, reco_pi0_subleading_photon_dir)
            reco_pi0_opening_angle = np.degrees(np.arccos(reco_pi0_cos_opening_angle))
            reco_pi0_mass = np.sqrt(2*reco_pi0_photons[0].calo_ke*reco_pi0_photons[1].calo_ke*(1-reco_pi0_cos_opening_angle))
            
            # reco dict corresponding to a CSV row
            reco_pi0_dict = {}
            reco_pi0_dict['reco_interaction_id'] = match[0].id
            reco_pi0_dict['reco_vertex_x'] = reco_inter.vertex[0]
            reco_pi0_dict['reco_vertex_y'] = reco_inter.vertex[1]
            reco_pi0_dict['reco_vertex_z'] = reco_inter.vertex[2]
            
            reco_pi0_dict['reco_leading_ph_start_point_x'] = reco_pi0_photons[0].start_point[0]
            reco_pi0_dict['reco_leading_ph_start_point_y'] = reco_pi0_photons[0].start_point[1]
            reco_pi0_dict['reco_leading_ph_start_point_z'] = reco_pi0_photons[0].start_point[2]
            reco_pi0_dict['reco_leading_ph_energy'] = reco_pi0_photons[0].calo_ke
            
            reco_pi0_dict['reco_subleading_ph_start_point_x'] = reco_pi0_photons[1].start_point[0]
            reco_pi0_dict['reco_subleading_ph_start_point_y'] = reco_pi0_photons[1].start_point[1]
            reco_pi0_dict['reco_subleading_ph_start_point_z'] = reco_pi0_photons[1].start_point[2]
            reco_pi0_dict['reco_subleading_ph_energy'] = reco_pi0_photons[1].calo_ke
            
            reco_pi0_dict['opening_angle'] = reco_pi0_opening_angle
            reco_pi0_dict['mass'] = reco_pi0_mass
            reco_pi0_dict['opening_angle'] = reco_pi0_opening_angle
            reco_pi0_dict['mass'] = reco_pi0_mass 

            ### To-do ##########################################################
            # Add truth info, so we can tell if we are selecting pi0s properly
            #
            #
            ####################################################################
            
            # Append row to CSV
            self.append('log', **reco_pi0_dict)
            
            
                
            
