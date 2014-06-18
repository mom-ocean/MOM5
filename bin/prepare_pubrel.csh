#!/bin/csh
#
#These are the notes on how to make the MOM4p1 public release.
#The notes are intended only for the GFDL people responsible for making the MOM4p1 public release.
#
set FMS_RELEASE=tikal
set MOM_RELEASE=tikal
set PUB_RELEASE=mom5_pubrel_December2013
mkdir src
cd src

#check outs
cvs co -r $MOM_RELEASE mom5
cvs co -r $FMS_RELEASE shared ocean_shared
cvs co -r $FMS_RELEASE ice_sis ice_param
cvs co -r $FMS_RELEASE land_lad land_null land_param land_lad2
cvs co -r $FMS_RELEASE atmos_param_am3
cvs co -r $FMS_RELEASE atmos_bgrid atmos_coupled atmos_ebm atmos_fv_dynamics
cvs co -r $FMS_RELEASE atmos_null atmos_shared atmos_spectral
cvs co -r $FMS_RELEASE coupler

cvs co -r $FMS_RELEASE postprocessing preprocessing  tools

#cut and nullify some stuff
\rm -rf atmos_param/qe_moist_convection atmos_param/two_stream_gray_rad atmos_param/shallow_physics 
\rm -rf atmos_param/lin_cloud_microphys/* atmos_param/clubb/*
cvs co -r nullify_rab_nnz atmos_param/lin_cloud_microphys/lin_cloud_microphys.F90 atmos_param/clubb/CLUBB_driver_SCM.F90 atmos_param/clubb/MG_microp_3D.F90

#The following file is not needed
\rm -rf mom5/drivers/coupler_types.*

\rm -rf tools/xmlDoc tools/fremetar tools/fbrowser


cvs up -r mom5_pubrel_dec2013_nnz mom5/doc
mv mom5/doc ../

#No pdf,ps,html
find . -name '*.pdf' -exec rm -f {} \;
find . -name '*.ps' -exec rm -f {} \;
find . -name '*.html' -exec rm -f {} \;

#bin/ and exp/

cvs co -r mom5_pubrel_dec2013_nnz mom5/utils/
mv mom5/utils/bin ../
mv mom5/utils/exp ../
\rm -rf mom5/utils

#Date tag (sticky) the whole thing
#cvs tag $(PUB_RELEASE)_nnz *

exit 
#The following cleanups is done after tagging.
#No www.gfdl.noaa.gov
foreach file ( `grep -l -r noaa.gov .` )
foreach? sed '/.*www.*.noaa.gov.*/d' -i $file
foreach? end

#No work emails
foreach file ( `grep -l -r @noaa.gov .`)
foreach? sed 's/@noaa.gov/@no.gov/g' -i $file
foreach? end

foreach file ( `grep -l -r @gfdl.noaa.gov .`)
foreach? sed 's/@gfdl.noaa.gov/@no.gov/g' -i $file
foreach? end

foreach file ( `grep -l -r EMAIL .`)
sed 's/EMAIL=.*.gov/EMAIL="GFDL.Climate.Model.Info@noaa.gov/g' -i  $file
end

foreach file ( `grep -l -r @no.gov .` )
foreach? sed 's/@.*no.gov//g' -i $file
foreach? end

exit

#update the htmls for some major files
#The following is too big to be included
cvs co -r xmlDoc_clean_sdu tools/xmlDoc


setenv XMLDOC_ROOT ${PWD}/tools/xmlDoc

foreach file ( coupler shared/diag_manager shared/mpp shared/time_manager shared/data_override preprocessing shared/field_manager mom5 )
${PWD}/tools/xmlDoc/bin/xmlDoc --dir $file
end

#GNUize all F90 and c files
#Some files are checked in without write permission
#chmod -R +w .
#/home/nnz/bin/GNULicense.pl -f --dir=. --recursive
#This will touch ALL code, so you have to test, test, test!
# cvs ci or not ci ,  I wouldn't do it GNU!!

#
#tar up
#
cd ../


tar cvf $PUB_RELEASE.tar $PUB_RELEASE/src
tar rvf $PUB_RELEASE.tar $PUB_RELEASE/bin
tar rvf $PUB_RELEASE.tar $PUB_RELEASE/exp
gzip $PUB_RELEASE.tar

#Documentation
#cvs co -r $MOM_RELEASE mom5/doc
#cd mom5/doc
##Work on the doc/README check it in and move it up
##Work on the doc/quickstart_guide.xml and MOM_practice.xml
##Generate the html and pdf from xml
#/home/nnz/bin/mkdocbk quickstart_guide.xml
#/home/nnz/bin/mkdocbk MOM_practice.xml
#Run Seth's tool to produce .html files for time_manager.html field_manager.html diag_manager.html mpp.html mpp_io.html coupler_main.html
#and move them in doc/  
#cvs ci README quickstart_guide.xml quickstart_guide.html MOM_practice.xml MOM_practice.html  time_manager.html field_manager.html diag_manager.html mpp.html mpp_io.html coupler_main.html
#mv README ../
