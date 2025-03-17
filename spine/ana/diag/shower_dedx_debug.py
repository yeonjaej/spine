'''Module to evaluate diagnostic metrics on showers.'''

import numpy as np

from spine.ana.base import AnaBase
from spine.utils.gnn.cluster import cluster_dedx_legacy, cluster_dedx_DBScan_PCA, cluster_dedx_dir
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist

__all__ = ['ShowerStartSingleParticle']

def cluster_dedx2_with_PCA_debug(voxels,
                 values,
                 start,
                 dedx_dist=3, cont_dist=5, detailed=True, num_intervals=10):
    # If max_dist is set, limit the set of voxels to those within a sphere of radius max_dist                                
    assert voxels.shape[1] == 3, (
            "The shape of the input is not compatible with voxel coordinates.")

    # If start point is not in voxels, assign the closest point within voxels
    # as the startpoint
    if not np.isclose(start, voxels, atol=1e-2).all(axis=1).any():
        dists = np.linalg.norm(voxels - start, axis=1)
        perm = np.argsort(dists)
        start = voxels[perm[0]]

    # distance from the startpoint
    dist_mat = cdist(start.reshape(1,-1), voxels).flatten()

    # legacy dedx
    #if dedx_dist > 0:
    voxels_dedx = voxels[dist_mat <= dedx_dist]
    #print("thr: ", max_dist, ", num vox: ", len(voxels))
    if len(voxels_dedx) < 2:
        return 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,[0., 0., 0.]
    values_dedx = values[dist_mat <= dedx_dist]
    dist_dedx = dist_mat[dist_mat <= dedx_dist]
    # Calculate sum of values                                                                                                                                                                                                         
    sum_dedx = np.sum(values_dedx)
    # Calculate max distance for dedx
    max_dist_dedx = np.max(dist_dedx)

    
    # continuity check 
    voxels_cont = voxels[dist_mat <= cont_dist]
    values_cont = values[dist_mat <= cont_dist]
    dist_cont = dist_mat[dist_mat <= cont_dist]
    # Perform DBSCAN clustering
    # parameters are not yet tuned
    eps = 0.59
    min_samples = 1
    dbscan = DBSCAN(eps, min_samples=min_samples)
    cluster_labels = dbscan.fit_predict(voxels_cont)
    #clusts, counts = np.unique(cluster_labels, return_counts=True)
    num_clust = max(0, max(cluster_labels)+1)

    # fining the dbscan cluster containing the startpoint
    start_clust = -1
    true_p_clust1 = -1
    true_p_clust2 = -1
    true_p_clust3 = -1
    for i in range(num_clust):
        if np.isclose(start, voxels_cont[cluster_labels==i], atol=1e-2).all(axis=1).any():
            start_clust = i
        if i==0:
            true_p_clust1 = len(voxels_cont[cluster_labels==i])
        if i==1:
            true_p_clust2 = len(voxels_cont[cluster_labels==i])
        if i==2:
            true_p_clust3 = len(voxels_cont[cluster_labels==i])
            
    voxels_clust = voxels_cont[cluster_labels==start_clust]
    values_clust = values_cont[cluster_labels==start_clust]
    dist_clust = dist_cont[cluster_labels==start_clust]
    
    voxels_clust = voxels_clust[dist_clust<=dedx_dist]
    values_clust = values_clust[dist_clust<=dedx_dist]
    dist_clust = dist_clust[dist_clust<=dedx_dist]

    if len(voxels_clust)<3:
        return sum_dedx, max_dist_dedx, 0., 0., 0., 0., 0., 0., 0., 0., len(voxels_clust), 0., 0., 0., [0., 0., 0.]
    dedx_clust_dist = np.max(dist_clust)
    # include dE from other clusters
    voxels_inc = voxels_cont[dist_cont<=dedx_clust_dist]
    values_inc = values_cont[dist_cont<=dedx_clust_dist]
    dist_inc = dist_cont[dist_cont<=dedx_clust_dist]

    # Perform PCA
    pca = PCA(n_components=3)
    pca.fit(voxels_clust)
    p_axis = pca.components_[0]
    p_fit = pca.explained_variance_ratio_[0]
    
    # Project voxels onto the principal axis
    p_voxels = np.dot(voxels_clust - np.mean(voxels_clust, axis=0), p_axis)

    min_proj = np.min(p_voxels)
    max_proj = np.max(p_voxels)
    mask = (p_voxels >= min_proj) & (p_voxels <= max_proj)

    # inclusive voxels
    p_voxels_inc = np.dot(voxels_inc - np.mean(voxels_clust, axis=0), p_axis)
    min_proj_inc = np.min(p_voxels_inc)
    max_proj_inc = np.max(p_voxels_inc)
    #mask_inc = (p_voxels_inc >= min_proj_inc) & (p_voxels_inc <= max_proj_inc)
    mask_inc = (p_voxels_inc >= min_proj) & (p_voxels_inc <= max_proj)
    #print("len: ", len(values_dedx), "proj len: ", len(values_dedx[mask]))
    p_sum = np.sum(values_clust[mask])

    voxels_de = voxels_inc[mask_inc]
    values_de = values_inc[mask_inc]
    p_sum_inc = np.sum(values_de)
    p_length = max_proj - min_proj
    p_length_inc = max_proj_inc - min_proj_inc

    voxels_sp = voxels_de - np.mean(voxels_clust, axis=0)
    p_voxels_sp = np.dot(voxels_sp, p_axis)
    vectors_to_axis = voxels_sp - np.outer(p_voxels_sp, p_axis)
    spread = np.linalg.norm(vectors_to_axis, axis=1)
    spread = sum(spread)/len(voxels_sp)
    #print("spread: ", spread)
    
    return sum_dedx, max_dist_dedx, p_sum, p_length, p_sum_inc, p_length, spread, p_fit, num_clust, start_clust, len(voxels_clust), true_p_clust1, true_p_clust2, true_p_clust3, p_axis

class ShowerStartSingleParticle(AnaBase):
    """This analysis script computes the dE/dx value within some distance
    from the start point of an EM shower object.

    This is a useful diagnostic tool to evaluate the calorimetric separability
    of different EM shower types (electron vs photon), which are expected to
    have different dE/dx patterns near their start point.
    """

    # Name of the analysis script (as specified in the configuration)
    name = 'shower_dedx_debug'
    
    def __init__(self, radius, is_photon=False, **kwargs):
        """Initialize the analysis script.

        Parameters
        ----------
        radius : Union[float, List[float]]
            Radius around the start point for which evaluate dE/dx
        **kwargs : dict, optional
            Additional arguments to pass to :class:`AnaBase`
        """
        # Initialize the parent class
        super().__init__('particle', 'both', **kwargs)

        # Store parameters
        self.radius = radius

        self.is_photon = is_photon

        # Initialize the CSV writer(s) you want
        for obj in self.obj_type:
            self.initialize_writer(obj)
        self.update_keys({'clust_label_adapt': True, 'meta': True, 'particles': True, 'clust_label_g4': True})
        self.units = 'cm'

    def process(self, data):
        """Evaluate shower start dE/dx for one entry.

        Parameters
        ----------
        data : dict
            Dictionary of data products
        """
        # Fetch the keys you want
        print("index: ", data['index'], ", len of truth: ", len(data['truth_particles']))
        if (len(data['truth_particles']) > 0):
            
            out_dict = {}
            out_dict['index'] = data['index']
            #print(data['index'])
            out_dict['file_index'] = data['file_index']
            out_dict['file_entry_index'] = data['file_entry_index']

            out_dict['num_true'] = len(data['truth_particles'])
            # True Showers

            true_shower = data['truth_particles'][0]

            out_dict['true_creation_process'] = true_shower.creation_process
                
            startpoint = true_shower.start_point
            
            match_id = -1
            if true_shower.match_overlaps is not None and len(true_shower.match_overlaps)>0:
                out_dict['true_match_overlap'] = true_shower.match_overlaps[0]
                match_id = true_shower.match_ids[0]
            # Reco Showers, reco points
            
            if len(data['reco_particles']) == 0:
                return
            if match_id == -1:
                return

            reco_shower = data['reco_particles'][match_id]            
            startpoint = reco_shower.start_point
            out_dict['reco_particle_id'] = reco_shower.id
            # out_dict['reco_particle_energy_init'] = reco_shower.energy_init
            out_dict['reco_particle_energy_deposit'] = reco_shower.calo_ke
            out_dict['reco_is_contained'] = reco_shower.is_contained
            out_dict['reco_match_overlap'] = -1
            if reco_shower.match_overlaps is not None and len(reco_shower.match_overlaps)>0:
                out_dict['reco_match_overlap'] = reco_shower.match_overlaps[0]
            sum_dedx, max_dist_dedx, p_sum, p_length, p_sum_inc, p_length, spread, p_fit, num_clust, start_clust, start_clust_size, p_clust1, p_clust2, p_clust3, p_axis = cluster_dedx2_with_PCA_debug(reco_shower.points, reco_shower.depositions, startpoint, dedx_dist=self.radius)
            reco_p_de_inc, reco_p_dx_inc, reco_spread, reco_p_fit, reco_num_clust, reco_start_clust_size, reco_p_axis, reco_clust_sizes  = cluster_dedx_DBScan_PCA(reco_shower.points, reco_shower.depositions, startpoint, dedx_dist=self.radius, detailed=True)
            #out_dict['reco_de'] = reco_de_1
            #out_dict['reco_dx'] = reco_dx_1
            #out_dict['reco_PCA_de'] = reco_p_de
            #out_dict['reco_PCA_dx'] = reco_p_dx
            out_dict['reco_PCA_de_inc'] = p_sum_inc
            out_dict['reco_PCA_dx_inc'] = p_length
            out_dict['reco_p_spread'] = reco_spread
            out_dict['reco_p_fit'] = reco_p_fit
            out_dict['reco_p_num_clust'] = reco_num_clust
            
            out_dict['reco_p_start_clust_size'] = start_clust_size
            out_dict['reco_p_clust_size1'] = p_clust1
            out_dict['reco_p_clust_size2'] = p_clust2
            out_dict['reco_p_clust_size3'] = p_clust3
            out_dict['reco_start_clust'] = start_clust
            
            self.append('particle', **out_dict)
