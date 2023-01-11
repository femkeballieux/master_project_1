#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Takes in all the seperate vlass catalog and turns them into one
"""
from astropy.table import vstack, Table
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord

#run on desktop
path = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images/catalog_outputs/'

def read_file(batch, filename):
    """
    Takes in one of the single catalog and puts it in the proper format
    """
    hdulist = Table.read(path +  batch + '/' + filename)
    return(hdulist)

#way to run over all batches
batches = os.listdir(path)
list_list = [] #Stores all the single tables

#Run over all batches
for i, batch in enumerate(batches):
    #run over all fils in a batch
    filenames = os.listdir(path + batch)
    for j, filename in enumerate(filenames):
        #Get the table
        table_list = read_file(batch, filename)
        
        #store which batch it is
        batch_info = np.full(len(table_list[0][:]), batch , dtype="S7" )
        table_list['batch']=batch_info
        
        #store the index, to rule out duplicates
        index_info = np.full(len(table_list[0][:]), filename[0:3], dtype="S5" )
        table_list['index'] = index_info
        
        #Store the observation date
        date_info = np.full(len(table_list[0][:]), filename[-20:-10] , dtype="S13" )
        table_list['date'] = date_info
        
        #Add it to the list keeping track of all lists
        list_list.append(table_list)

#This is the final catalog
catalog = vstack(list_list)


catalog.write('vlass_catalog.fits', overwrite = True)
    

#TODO: do the indexing> right now it does not take into account the different lengths of the indexes, maybe something with 
#all the numbers before _ or something

