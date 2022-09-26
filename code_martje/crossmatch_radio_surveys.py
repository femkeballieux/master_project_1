# !/bin/python
# Crossmatching script for LoTSS sub-sample and other radio surveys.

import os

#When working on vdesk
#os.system('topcat -stilts tmatchn matcher=sky multimode=pairs nin=7 params=15 \
#    in1=/net/vdesk/data2/bach1/ballieux/master_project_1/data/unresolved_isolated_S_source.fits  values1="RAJ2000 DEJ2000" \
#    in2=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/NVSS.fits values2="RAJ2000 DEJ2000" \
#    in3=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/TGSS.fits values3="RAJ2000 DEJ2000" \
#    in4=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/VLSSr.fits values4="RAJ2000 DEJ2000" \
#    in5=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/LoLSS.fits values5="RA DEC"\
#    in6=/net/vdesk/data2/bach1/ballieux/master_project_1/data/first_14dec17.fits values6="RA DEC"\
#    in7=/net/vdesk/data2/bach1/ballieux/master_project_1/data/surveys/inband_spec_LoTSS.fits values7="RA DEC" \
#    out=/net/vdesk/data2/bach1/ballieux/master_project_1/data/crossmatch_NVSS_TGSS_VLSSr_LoLSS_FIRST_inband.fits')

#First was removed?

#When working on laptop
os.system('C:/Users/Femke/Documents/GitHub/master_project_1/code_martje/topcat-full.jar -stilts \
          tmatchn matcher=sky multimode=pairs nin=7 params=15 \
    in1=C:/Users/Femke/Documents/GitHub/master_project_1/data/unresolved_isolated_S_source.fits  values1="RA DEC" \
    in2=C:/Users/Femke/Documents/GitHub/master_project_1/data/surveys/NVSS.fits values2="RAJ2000 DEJ2000" \
    in3=C:/Users/Femke/Documents/GitHub/master_project_1/data/surveys/TGSS.fits values3="RAJ2000 DEJ2000" \
    in4=C:/Users/Femke/Documents/GitHub/master_project_1/data/surveys/VLSSr.fits values4="RAJ2000 DEJ2000" \
    in5=C:/Users/Femke/Documents/GitHub/master_project_1/data/surveys/LoLSS.fits values5="RA DEC"\
    in6=C:/Users/Femke/Documents/GitHub/master_project_1/data/surveys/inband_spec_LoTSS.fits values7="RA DEC" \
    out=C:/Users/Femke/Documents/GitHub/master_project_1/data/crossmatch_NVSS_TGSS_VLSSr_LoLSS_inband.fits')