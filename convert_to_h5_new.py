#!/usr/bin/python
from __future__ import print_function, division
import os
import argparse
import uproot
import numpy as np
import h5py
import awkward1 as ak
#import akward as ak

def to_np_array(ak_array, maxN=100, num_entries = 576, pad=0):
    '''convert awkward array to regular numpy array'''
    print(ak_array)
    output0 = ak.pad_none(ak_array,num_entries,clip=True,axis=0)
    print("pad 0--------")
    #print(output0)
    output1 = ak.pad_none(output0,maxN,clip=True,axis=1)
    print("pad 1--------")
    #print(output1)
    output2 = ak.fill_none(ak.to_list(output1),pad)
    print("pad 2--------")
    #print(output2)
    return ak.to_numpy(output2)
#    return ak.fill_none(ak.pad_none(ak_array,maxN,clip=True,axis=-1),pad).to_numpy()

def store_objects(arrays, nentries, nobj=10, obj='jet'):
    '''store objects in zero-padded numpy arrays'''
    l1Obj_cyl = np.zeros((nentries,nobj,3))
    l1Obj_cart = np.zeros((nentries,nobj,3))
    pt = to_np_array(arrays['{}Et'.format(obj)],maxN=nobj, num_entries=nentries)
    eta = to_np_array(arrays['{}Eta'.format(obj)],maxN=nobj,num_entries=nentries)
    phi = to_np_array(arrays['{}Phi'.format(obj)],maxN=nobj,num_entries=nentries)
    #print(pt)
    #pt.flatten()
    #print("Flattening--------------------------------------------")
    #print(pt)
    l1Obj_cyl[:,:,0] = pt
    l1Obj_cyl[:,:,1] = eta
    l1Obj_cyl[:,:,2] = phi
    l1Obj_cart[:,:,0] = pt*np.cos(phi)
    l1Obj_cart[:,:,1] = pt*np.sin(phi)
    l1Obj_cart[:,:,2] = pt*np.sinh(eta)

    # now sort in descending pT order if needed
    sort_indices = np.argsort(-pt,axis=1)
    check_indices = np.tile(np.arange(0,nobj),(pt.shape[0],1))
    if not np.allclose(sort_indices,check_indices):
        l1Obj_cyl[:,:,0] = np.take_along_axis(l1Obj_cyl[:,:,0], sort_indices, axis=1)
        l1Obj_cyl[:,:,1] = np.take_along_axis(l1Obj_cyl[:,:,1], sort_indices, axis=1)
        l1Obj_cyl[:,:,2] = np.take_along_axis(l1Obj_cyl[:,:,2], sort_indices, axis=1)
        l1Obj_cart[:,:,0] = np.take_along_axis(l1Obj_cart[:,:,0], sort_indices, axis=1)
        l1Obj_cart[:,:,1] = np.take_along_axis(l1Obj_cart[:,:,1], sort_indices, axis=1)
        l1Obj_cart[:,:,2] = np.take_along_axis(l1Obj_cart[:,:,2], sort_indices, axis=1)
    return l1Obj_cyl, l1Obj_cart

def convert_to_h5(input_file, output_file, tree_name):
    inFile = uproot.open(input_file)
    #l1NtupleProduce/Stage3Regions/RegionTree
    l1Tree = inFile[tree_name]
    #nentries = l1Tree.num_entries
    ##fixme
    #nentries = 576


    # index of sum to save (MET)
    metindex = 2

    # save up to n objects (jets, muons, electrons)
    njets = 10
    nmuons = 4
    nelectrons = 4
    nregions = 576

    cylNames = [b'pT', b'eta', b'phi']
    cartNames = [b'px', b'py', b'pz']

    # variables to retrieve

    varList = ['vRegionEt','vRegionEta','vRegionPhi','vRegionEG','vRegionTau']

    #varList = ['nSums','sumEt','sumPhi',
    #           'jetEt', 'jetEta', 'jetPhi',
    #           'muonEt', 'muonEta', 'muonPhi',
    #           'egEt', 'egEta', 'egPhi']
    
    # get awkward arrays
    arrays = l1Tree.arrays(varList)

    # sums: store the following                                                                              
    # kTotalEt, kTotalEtEm, kTotalHt, kMissingEt, kMissingHt,                                                
    # with type 0, 16, 1, 2, 3
    #l1sum_cyl = np.zeros((nentries,3))
    #l1sum_cart = np.zeros((nentries,3))
    #sumEt = to_np_array(arrays['sumEt'],maxN=arrays['nSums'][0])
    #sumPhi = to_np_array(arrays['sumPhi'],maxN=arrays['nSums'][0])
    #l1sum_cyl[:,0] = sumEt[:,metindex] # MET_pt
    #l1sum_cyl[:,2] = sumPhi[:,metindex] # MET_phi
    #l1sum_cart[:,0] = sumEt[:,metindex]*np.cos(sumPhi[:,metindex]) # MET_px
    #l1sum_cart[:,1] = sumEt[:,metindex]*np.sin(sumPhi[:,metindex]) # MET_py

    # store objects: jets, muons, electrons
    ## check me
    nentries = len(arrays['vRegionEt'])
    print("nentries: " + str(nentries))

    l1Region_cyl, l1Region_cart = store_objects(arrays, nentries, nobj=nregions, obj='vRegion')
    
    print("cyl data size: " + str(l1Region_cyl.shape))
    #l1Jet_cyl, l1Jet_cart = store_objects(arrays, nentries, nobj=njets, obj='jet')
    #l1mu_cyl, l1mu_cart = store_objects(arrays, nentries, nobj=nmuons, obj='muon')
    #l1ele_cyl, l1ele_cart = store_objects(arrays, nentries, nobj=nelectrons, obj='eg')

    with h5py.File(output_file, 'w') as outFile:
        #outFile.create_dataset('FeatureNames_cyl', data=cylNames, compression='gzip')
        #outFile.create_dataset('FeatureNames_cart', data=cartNames, compression='gzip')
        outFile.create_dataset('l1Region_cyl', data=l1Region_cyl, compression='gzip')
        outFile.create_dataset('l1Region_cart', data=l1Region_cart, compression='gzip')
        #outFile.create_dataset('l1Jet_cyl', data=l1Jet_cyl, compression='gzip')
        #outFile.create_dataset('l1Jet_cart', data=l1Jet_cart, compression='gzip')
        #outFile.create_dataset('l1Muon_cyl', data=l1mu_cyl, compression='gzip')
        #outFile.create_dataset('l1Muon_cart', data=l1mu_cart, compression='gzip')
        #outFile.create_dataset('l1Ele_cyl', data=l1ele_cyl, compression='gzip')
        #outFile.create_dataset('l1Ele_cart', data=l1ele_cart, compression='gzip')
        #outFile.create_dataset('l1Sum_cyl', data=l1sum_cyl, compression='gzip')
        #outFile.create_dataset('l1Sum_cart', data=l1sum_cart, compression='gzip')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, required=True)
    parser.add_argument('--output-file', type=str, required=True)
    parser.add_argument('--tree-name', type=str, default='l1UpgradeEmuTree/L1UpgradeTree')
    args = parser.parse_args()
    convert_to_h5(**vars(args))
