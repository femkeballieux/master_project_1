"""Code to run AEGEAN"""
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

"""
This code is run in virtualenv

Ga naar VLASS_Aegean  mapje  
bash
cd /venv/bin
source activate

dan cd ..
cd ..

hiermee ben je terug in vorige mapje maar wel in venv. Vervolgens kan je gewoon Aegean_script runnen

voorbeeld van werkende code in command line is 
aegean --autoload --slice 0 --seedclip 3 --table 65_out.fits /net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/Tests/seedclip_4/065-VLASS_2020-06-26_J131747.37+602934.5_s1.0arcmin_e21_mosaicked.fits

"""

def classifyVlassFits(fitsName):
  """
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



def get_download(ra, dec, radius = 0.25 * u.arcmin, index=''):
    """
    Downloads the fitsfile for some coordinates, for a given radius
    """
    start=time.time()
    save_to_path ='/net/vdesk/data2/bach1/ballieux/master_project_1/VLASS_Aegean/VLASS_images_all/'
    coordinate = SkyCoord(ra, dec, frame='icrs', unit='deg')
    (CadcUrl,vlassTypes, fitsnames) = getCadcVlassUrl(coordinate, radius, epoch=["2.1", "2.2"], imageType='ql.tt0')
    print(coordinate, "\n", vlassTypes,"\n",CadcUrl,"\n",fitsnames, "\n" )
    for i, fitsname in enumerate(fitsnames): #possible that multiple fitsfiles get returned
        urllib.request.urlretrieve(CadcUrl[i], filename= save_to_path + index + fitsname)
    end=time.time()
    print('elapsed time for this source is  {:.5} s'.format(end-start))
#get_download(161.804, 47.0591)

path_to_mastersample= '/net/vdesk/data2/bach1/ballieux/master_project_1/data/'
mastersample = 'master_LoLSS_with_inband.fits'

hdulist = fits.open(path_to_mastersample + mastersample)
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
#print(orig_cols)
hdulist.close()

RA = tbdata['RA']

Dec =tbdata['Dec']

start_full = time.time()
for i, coords in enumerate(RA[:3]):
    get_download(RA[i], Dec[i], index=str(i)+'_')
    
end_full=time.time()
print('total elapsed time is {:.5} s'.format(end_full-start_full))

"""
The following code was written for doing only 767 PS sources, so commented it out. Will first try to obtain the vlass data in another way. 
Then once that works, in this same script run BANE and aegean. If that works, can we do in the same script the making of the catalogue?
Working with functions, need something to capture the time,exceptions?
"""
# #We split the process up into different patches
# path = 'VLASS_images/'  
# batch = '402_501/'

# #list all the filenames, unfortunatwly not in order, but the first numbers correspond to the index in the coord file
# filenames = os.listdir(path + batch) #only works once since after this we get 3 files for every folder. So need to duplicate the data maybe

# #The only way we can get the coordinates is where we inputted them
# coordinate_file = genfromtxt(path + 'PS_coords_402_501.csv', delimiter=',')

# #get the indexes
# index_array=np.zeros(len(filenames))
# for i, name in enumerate(filenames): #go trough all the filenames in a particular folder
#     index_array[i]=filenames[i][0:3] #This one contains the first 3 characters of every filename, not in order, contains duplicates
    
# #now we want to mask those that are duplicate, since these have measurements at different times
# m = np.zeros_like(index_array, dtype=bool) 
# m[np.unique(index_array, return_index=True)[1]] = True
# print('Store this number to check the duplicates',index_array[~m]) #These are the numbers of the duplicates

# #run bane and aegean
# for i, index in enumerate(filenames):
#     print("")
#     print('RUNNING BANE for', i, '/', len(filenames), 'which is', int(index_array[i]), 'in Cirada')
#     print("")
#     os.system('/home/ballieux/.local/bin/BANE --cores 1 ' + path + batch + filenames[i])
    
#     #load in the coordinates
#     coords =coordinate_file[int(index_array[i])]
#     RA= coords[0]
#     Dec= coords[1]
    
#     print("")
#     print('RUNNING aegean for', i, '/', len(filenames))
#     print("")
#     os.system('/home/ballieux/.local/bin/aegean --autoload --slice 0 --table ' + path + 'catalog_outputs/' + batch + str(int(index_array[i])) \
#               + '_'+ str(RA)[0:7] + '+' + str(Dec)[0:7]+ '_' + str(filenames[i])[10:20] + '.fits ' + path + batch + filenames[i] )

