
"""Code downloads vlass cutouts from cadc server, and runs BANE and aegean on it"""

import numpy as np
import os
from numpy import genfromtxt
from astropy.wcs.wcs import NoConvergence
from astroquery.cadc import Cadc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import table  
import urllib
import time
from astropy.io import fits
from astropy.table import vstack, Table
import matplotlib.pyplot as plt
from astropy.io import fits
import multiprocessing 
from functools import partial



"""
This code is run in my virtualenv

Ga naar VLASS_Aegean mapje  
bash
cd /venv/bin
source activate

dan cd ..
cd ..

hiermee ben je terug in vorige mapje maar wel in venv. Vervolgens kan je gewoon Aegean_script runnen
"""

def classifyVlassFits(fitsName):
  """
  Written by G. R. Sivakoff
  Given a VLASS FITS file following standard naming conventions, determine the type of image the FITS file is.
  All VLASS FITS files except spectral index (alpha) maps should follow *pbcor.AAA.subim.fits, 
  where AAA can be a variable length string that indicates the file type. VLASS FITS files that are spectral index
  (alpha) maps should follow *alpha.AAA.subim.fits.
  There does not seem to be a simple CAOM query that can do this, which is why we are using the name.
  """
  if ("pbcor" in fitsName): #Most files
    suffix = fitsName.split('pbcor.')
  elif ("alpha" in fitsName): #alpha files, we don' t really use these
    suffix = fitsName.split('alpha.')
    suffix[1] = "alpha."+suffix[1]
  if len(suffix) > 1:
    suffix = suffix[1]
  else:
    return(None)
  if len(suffix)>13:
    vlassType=suffix[:-11]
  else:
    vlassType=None
  return(vlassType)

def getVlassFitsNames(cadc,cadcResults):
  """
  Written by G. R. Sivakoff
   Given a CADC astroquery connection and the return of an already
   executed base query, get the names and set the types of all FITS files
   and classify their type within VLASS
   1. Search all auxillary files for an observation for fits files
   2. Classify FITS file type
  """
  fitsNames = np.array([])
  vlassTypes = np.array([])
  urls = cadc.get_data_urls(cadcResults,include_auxiliaries=True)
  for url in np.sort(urls):
    if (url[-4:] =='fits'):
      tempFitsName = url.rsplit('VLASS/')[1]
      fitsNames = np.append(fitsNames,tempFitsName)
      vlassTypes = np.append(vlassTypes,classifyVlassFits(tempFitsName))
  return (fitsNames, vlassTypes)

def getCoordVlass(targetCoordinate,targetRadius,
                  targetEpochs,targetTypes,
                  cadc):
  """
  Written by G. R. Sivakoff
   Given a target coordinate, radius, Epochs, and image types, this
   uses an active CADC astroquery connection that has already been
   verified to have results to return the URL for VLASS cutout download 
   and the image type among VLASS options.
  """
  fitsNames         = np.array([])
  vlassTypes        = np.array([])
  CadcUrl           = np.array([])
  completedSubtiles = np.array([])
  
  """
   Note that the following warnings show up everytime this is run
   WARNING: UnitsWarning: The unit 'pix' has been deprecated in the 
     VOUnit standard. [astropy.units.format.utils]
   WARNING:astroquery:UnitsWarning: The unit 'pix' has been deprecated in
     the VOUnit standard.
   /usr/local/lib/python3.7/dist-packages/astropy/table/table.py:3401: 
     FutureWarning: elementwise == comparison failed and returning scalar
     instead; this will raise an error or perform elementwise comparison
     in
     the future.
       result = (self.as_array().data == other) & 
          (self.mask == false_mask)
   /usr/local/lib/python3.7/dist-packages/astropy/table/table.py:3409:
     FutureWarning: elementwise == comparison failed and returning scalar
     instead; this will raise an error or perform elementwise comparison
     in
  """

  results = cadc.query_region(targetCoordinate, targetRadius,
                              collection='VLASS')

  if len(results['observationURI']) == 0: #This part I changed, checks if anything is found
    print("No Valid Coordinate")
    return(None)
  else:
    # Reject QA failed observations
    passedResults =   results[results['requirements_flag']!='fail']
    # Sort by obstime
    sortedSubtiles = (passedResults.group_by('time_bounds_lower'))
    for subtile in sortedSubtiles:
      # This needs to be done becausethe same observationURI can
      # have multiple results. In particular, one for Quick Look
      # (currently calibrationLevel==2) and one for Single Epoch
      # (currently calibrationLevel==3). This is not ideal, but, ...
      if subtile['observationURI'] not in completedSubtiles:
        completedSubtiles = np.append(completedSubtiles,
                                      subtile['observationURI'])
        # Gets the epoch of the current subtile result
        observedEpoch = subtile['proposal_id'][5:]
        # Gets the observationURI of the current subtile result        
        observationURI = subtile['observationURI']
        # Only return urls for targetted Epochs
        if(observedEpoch in targetEpochs):
          # Since I think QL should always be there [we should check this
          # with NRAO], get all FITS files associated with an 
          # observationURI by first only considering the QL subtile 
          # results. This is done so one does not get the same files
          # twice, once for the astroquery result for QL and one for the
          # astroquery result for SE.
          resultSelection = np.logical_and(
              results['observationURI']==observationURI,
              results['calibrationLevel']==2)
          subtileResult = results[resultSelection]
          # Gets the VLASS FITS names for a QL observationURI
          (tempFitsNames,tempVlassTypes) = getVlassFitsNames(cadc,
                                                            subtileResult)
          # This is not an ideal for loop. It should probably be converted
          # to np.where
          # Ideally CADC should only have one set of QL and one set of SE
          # images per ObservationURI. There is no real error checking
          # being done about this...
          for ii in np.arange(len(tempVlassTypes)):
            if(tempVlassTypes[ii]=="tt0"):
              if("ql" in tempFitsNames[ii]):
                baseUrl = cadc.get_image_list(subtileResult,
                                              targetCoordinate,
                                              targetRadius)
                baseUrl = baseUrl[0]
                qlName = tempFitsNames[ii]
          # Now that we have the Quick Look default result, consider
          # all results with the same observationURI
          resultSelection = results['observationURI']==observationURI,
          subtileResult = results[resultSelection]
          # Gets the VLASS FITS names for a QL or SE observationURI
          (tempFitsNames,tempVlassTypes) = getVlassFitsNames(cadc,
                                                            subtileResult)
          for ii in np.arange(len(tempVlassTypes)):
            url = np.char.replace(baseUrl,qlName,tempFitsNames[ii])
            # This probably should use some better error catching.
            tempVlassType = "NOSUCHTYPE"
            # Prepend SE versus QL
            if("se" in tempFitsNames[ii]):
                tempVlassType = "se."+tempVlassTypes[ii]
            elif ("ql" in tempFitsNames[ii]):
                tempVlassType = "ql."+tempVlassTypes[ii]
            # Append values of CADC URL, FITS file names (not currently 
            #used),
            # and VLass Image Types
            if (tempVlassType in targetTypes):
              CadcUrl    = np.append(CadcUrl,url)    
              fitsNames  = np.append(fitsNames,tempFitsNames[ii])
              vlassTypes = np.append(vlassTypes,tempVlassType)       
    return(CadcUrl,vlassTypes, fitsNames)

def getCadcVlassUrl(coordinate, radius, 
                    epoch=np.array(["1.1","1.2",
                                    "2.1","2.2",
                                    "3.1","3.2"]),
                    imageType = np.array(['ql.tt0',
                                          'ql.tt0.rms',
                                          'se.tt0',
                                          'se.tt0.rms',
                                          'se.tt1',
                                          'se.tt1.rms',
                                          'se.alpha.error',
                                          'se.alpha'])):
  """
  Written by G. R. Sivakoff
   Given a target coordinate, radius, Epochs, and image types, this
   uses an active CADC astroquery connection to first verify that there
   will be results, and then to return the URL for VLASS cutout download 
   and the image type among VLASS options.
  """
  cadc = Cadc()
  inVlassCheck = checkCoordVlass(coordinate,epoch,cadc)
  if True in inVlassCheck:
    return(getCoordVlass(coordinate,radius,epoch,imageType,cadc))
  return(None,None)

def checkCoordVlass(targetCoordinate,targetEpochs,cadc):
  """
  Written by G. R. Sivakoff
   Given a target coordinate and Epochs, this
   uses an active CADC astroquery connection to verify that there
   will be results.
  """
  checkPassed =  np.array([], dtype='bool')

  # 2 arcsecond availability for quick check
  existenceRadius = 2/3600.*u.deg
  """
   Note that the following warnings show up everytime this is run
   WARNING: UnitsWarning: The unit 'pix' has been deprecated in the 
     VOUnit standard. [astropy.units.format.utils]
   WARNING:astroquery:UnitsWarning: The unit 'pix' has been deprecated in
     the VOUnit standard.
   /usr/local/lib/python3.7/dist-packages/astropy/table/table.py:3401: 
     FutureWarning: elementwise == comparison failed and returning scalar
     instead; this will raise an error or perform elementwise comparison
     in
     the future.
       result = (self.as_array().data == other) & 
          (self.mask == false_mask)
   /usr/local/lib/python3.7/dist-packages/astropy/table/table.py:3409:
     FutureWarning: elementwise == comparison failed and returning scalar
     instead; this will raise an error or perform elementwise comparison
     in
     the future. result = self.as_array() == other
  """
  results = cadc.query_region(targetCoordinate, existenceRadius,
                              collection='VLASS', verbose=True)
  if len(results['observationURI']) == 0: #Changed from original, checks if anything is found
    print("No Valid Coordinate")
    checkPassed = np.append(checkPassed,False)
  else:
    # Sort by obstime
    sortedSubtiles = (results.group_by('time_bounds_lower'))
    for subtile in sortedSubtiles:
      observedEpoch = subtile['proposal_id'][5:]
      if(observedEpoch in targetEpochs):
        checkPassed = np.append(checkPassed,True)
      else:
        checkPassed = np.append(checkPassed,False)
    return(checkPassed)




def get_url(coordinate, radius = 0.25 * u.arcmin):
    """
    Written by Femke Ballieux
    Retrieves the urls for the fitsfile for some coordinates, for a given radius.
    For epoch 2.1 and 2.2 only, only quick look images
    Index is included to keep track of different files, as well as to prevent files
    from overwriting eachother
    """
    start=time.time() 


    #Code below retrieves the urls, image type and filenames. Can be more than 1 image per query
    (CadcUrl, vlassTypes, fitsnames) = getCadcVlassUrl(coordinate, radius, epoch=["2.1", "2.2"], imageType='ql.tt0') 
    print(coordinate, "\n", vlassTypes,"\n",CadcUrl,"\n",fitsnames, "\n" )
    end=time.time()
    print("Retrieving urls in {:.5} seconds".format(end-start))
    return vlassTypes, CadcUrl, fitsnames

def get_download(coordinate, vlassTypes, CadcUrl, fitsnames, index=''):
    """
    Written by Femke Ballieux
    downloads the urls for a specific coordinate, using the output from the get-url function
    """
    start=time.time() 
    save_to_path ='/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images_all/' #Where the image will be stored
    for i, fitsname in enumerate(fitsnames): #possible that multiple fitsfiles get returned, run over all
        urllib.request.urlretrieve(CadcUrl[i], filename= save_to_path + index + fitsname) #retrieve the file
        
    end=time.time()
    print('Downloading urls in  {:.5} seconds'.format(end-start))
    return fitsnames
    
def run_BANE(filename, index):
    # begin_number = 0
    # end_number = 10000
    
    """
    Written by Femke Ballieux
    Runs BANE for a downloaded image. The index is an argument, which gets included 
    to deal with the proper filenames
    """
    path_to_image = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images_all/'
    # print('RUNNING BANE for index ', begin_number + index, '/', end_number-1)
    print('RUNNING BANE for index ',  index)
    print("")
    os.system('BANE --cores 1 ' + path_to_image + str(index) + '_' + filename)

def run_Aegean(filename, index, RA='_', Dec='_'):
    # begin_number = 0
    # end_number = 10000
    
    """
    Runs Aegean for the vlass images where BANE had already run on. Index and RA, Dec are included 
    to get the proper output filenames. 
    --autoload
    --seedclip 4 since the positions are known.
    --slice 0 to prevent any false cubes to interfere with the results
    """
    path_out = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_all_output/'
    path_in = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images_all/'
    print("")
    print('RUNNING aegean for index ', index)
    print("")

    #Output filename is index_RA+Dec_originalfilename.fits
    os.system('aegean --autoload --slice 0 --cores 1 --seedclip 4 --table ' + path_out + str(index) + '_' \
               + filename + '.fits ' \
            + path_in + str(index) + '_' + filename )
    # os.system('aegean --autoload --slice 0 --seedclip 4 --table ' + path_out + str(index) + '_' \
    #           + str(RA) + '+' + str(Dec) + '_' + filename + '.fits ' \
    #         + path_in + str(index) + '_' + filename )


# #Below I will do some testing/messing around. Delete below here
# path_to_mastersample= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/' 
# mastersample = 'master_LoLSS_with_inband.fits' #Mastersample is taken now

# #Read in our mastersample for the coordinates
# hdulist = fits.open(path_to_mastersample + mastersample)
# tbdata = hdulist[1].data
# orig_cols = hdulist[1].columns
# hdulist.close()

# #Coordinates
# RA = tbdata['RA']
# Dec =tbdata['Dec']

# test=5
# #result_table = Simbad.query_region(SkyCoord(ra=[10, 11], dec=[10, 11], unit=(u.deg, u.deg), frame='fk5'), radius=0.1 * u.deg)
# #result_table = Cadc.query_region(SkyCoord(ra=[10, 11], dec=[10, 11], unit=(u.deg, u.deg), frame='fk5'), radius=0.1 * u.deg)
# #print(result_table, 'hoi')

# coordinates = SkyCoord(ra=RA[:test], dec=Dec[:test], frame='icrs', unit='deg')
# print(coordinates)
# urls=getCadcVlassUrl(coordinates, radius=0.25*u.arcmin)
# #urls=Cadc.query_region(coordinates, radius=0.1*u.deg, collection='VLASS')
# print(urls)
# #What I learned so far, the extended version is better since it deals with filenames, epochs etc. Do need to figure out how to do batches. 
# # #Wat ik hierboven heb is om te testen zodat ik niet mn hele code verneuk 
# #Hij pakt niet meerdere coordinates in 1x. Misschien rondspelen met hoe je coordinates defined, het zou toch moeten kunnen?? 
# #Delete above here



def obtain_and_clean_images(coordinate, index = ''):
  """
  Obtains vlass url and downloads the images for a single coordinate, then cleans the downloaded images
  """
  start_full = time.time()

  #Get the url and the download
  try:
      vlassTypes, CadcUrl, fitsnames = get_url(coordinate)
      get_download(coordinate, vlassTypes, CadcUrl, fitsnames, str(index)+'_')
  except:
      print('DOWNLOADING WENT WRONG FOR', index)
          

  end=time.time()
  print('elapsed time for retrieving the url and downloading the image is {:.5} s'.format(end-start_full))
  
  
  for a, fitsname in enumerate(fitsnames): #Runs over multiple images for a single coordinate pair = index
    start_BANE=time.time()
    try:
        run_BANE(fitsnames[a],index) #Gets BANE for each image
    except:
        print('COULD NOT RUN BANE FOR', index, a)

    end_BANE=time.time()
    print('Doing BANE in {:.5} seconds'.format(end_BANE-start_BANE))

    start_Aegean=time.time()
    try:
        run_Aegean(fitsname, index) #Gets Aegean for each image
    except:
        print('COULD NOT RUN AEGEAN FOR', index, a)
    end_Aegean=time.time()
    print('Doing Aegean in {:.5} seconds'.format(end_Aegean-start_Aegean))
    print("")
    
  end_full=time.time()
  print('total elapsed time for this source is {:.5} s'.format(end_full-start_full))

    
def multiprocessing_func(coordinates, i):
    obtain_and_clean_images(coordinates[i], index=i)


    
"""
Below is the code that turns the results into a catalogue
"""

def read_file(filename):
    """
    Takes in one of the single catalog and puts it in the proper format
    """
    path_out = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_all_output/'
    hdulist = Table.read(path_out + '/' + filename)
    return(hdulist)

def make_catalog():
  """
  Reads in the output files from aegean and turns them into a single catalog
  """
  path_out = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_all_output/'
  list_list=[]#A list that keeps track of all the different tables

  #run over all files that have been outputted
  filenames = os.listdir(path_out)
  for j, filename in enumerate(filenames): #Run over all files
    print(j)
    #Get the table
    table_list = read_file(filename)

    #store the index, to rule out duplicates
    inx=filename.split('_')
    index_info = np.full(len(table_list[0][:]), inx[0], dtype="S5" ) #make an array with the index
    table_list['index'] = index_info
          
    #Add it to the list keeping track of all lists
    list_list.append(table_list)
  print('Doing the vstack')
  #This is the final catalog
  catalog = vstack(list_list)
  print('writing the catalog')
  catalog.write('vlass_catalog_all.fits', overwrite = True)

  """
  This function works fine for a few files, but it not good with memory for many files.
  what to do:
      cd to the directory of your aegean catalogs
      ls > filenames.txt
      move this file one directory up
      check if filenames.txt does not contain itself
      
      java -jar /net/vdesk/data2/bach1/ballieux/master_project_1/topcat-full.jar -stilts tcat in=@/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/filenames.txt out=/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/vlass_catalog.fits lazy=true seqcol=input_table uloccol=filename

      This function does not give an output, does it a few minutes
  """


if __name__ == '__main__':
 	# set the start method
    multiprocessing.set_start_method('spawn')
    # begin_number = 0
    # end_number = 10000
    
    """
    Below, we import the coordinates and turn it into a coordinate object. These coordinates
    are looped over for the downloading of the images
    """
    path_to_mastersample= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/' 
    path_out = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_all_output/'
    path_image = '/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images_all/'
    #mastersample = 'master_LoLSS_with_inband.fits' #Mastersample is taken now
    LoTSS_sample = 'crossmatch_NVSS_LoTSS.fits'

    #Read in our mastersample for the coordinates
    hdulist = fits.open(path_to_mastersample + LoTSS_sample)
    tbdata = hdulist[1].data
    orig_cols = hdulist[1].columns
    hdulist.close()

    #Coordinates
    RA = tbdata['RA']
    Dec = tbdata['DEC'] #Pay attention to proper column names

    coordinates = SkyCoord(ra=RA, dec=Dec, frame='icrs', unit='deg')
    
    #check which files are already cleaned by Aegean
    print('Checking which files already have an Aegean output')
    filenames = os.listdir(path_out)
    index_done=[int(filenames[i].split('_')[0]) for i in range(len(filenames))]
  #   index_done_array=np.array(index_done)

    mask_notdone = np.ones((len(coordinates)), dtype=bool)
    mask_notdone[index_done] = False
    
    print('Checking which files already have a cutout')
    filenames_image = os.listdir(path_image)
    index_image=[int(filenames_image[i].split('_')[0]) for i in range(len(filenames_image))]
#    index_image_array=np.array(index_image)
    mask_notdone[index_image] = False

    indexes=np.array(range(len(coordinates)))
    
    starttime = time.time()
    pool = multiprocessing.Pool()
    func = partial(multiprocessing_func, coordinates)
    pool.map(func, indexes[mask_notdone])
    pool.close()
    pool.join()
    print('That took {} seconds'.format(time.time() - starttime))
    
   #print('making catalog')
    # make_catalog() #Makes catalog of all sources found