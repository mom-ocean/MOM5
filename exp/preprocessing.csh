#!/bin/csh -f
#
#Modify the namelists in input.nml to comply with the changes in this release
#
  sed "s/barotropic_time_stepping_mom4p0/barotropic_time_stepping_A/g" -i input.nml
  sed "s/barotropic_time_stepping_mom4p1/barotropic_time_stepping_B/g" -i input.nml
  sed "/.*robert_asselin_bt.*/d" -i input.nml
  sed "/.*barotropic_leap_frog.*/d" -i input.nml
  sed "/.*barotropic_pred_corr.*/d" -i input.nml
  sed "s/linear_eos/eos_linear/g" -i input.nml
  sed "s/alpha_eos_linear/alpha_linear_eos/g" -i input.nml
  sed "s/beta_eos_linear/beta_linear_eos/g" -i input.nml
  sed "s/eos_linear=.false./eos_linear=.false. \n eos_preteos10=.true./g" -i input.nml
  sed "/.*freezing_temp_accurate.*/d" -i input.nml
#
#Correct field_table enteries for ppm_hlimiter=3 and ppm_vlimiter=3 to be 2 instead of 3
#
  sed "s/ppm_hlimiter=3/ppm_hlimiter=2/g" -i field_table
  sed "s/ppm_vlimiter=3/ppm_vlimiter=2/g" -i field_table
# 
#Avoid am3 physics
#
#  sed "s/do_lsc          =.true./do_lsc          =.false./g" -i input.nml
#  sed "s/do_strat        =.false./do_strat        =.true./g" -i input.nml
#  sed "s/do_diffusivity   = .true./do_diffusivity   = .false./g" -i input.nml
#
#MOM5 required changes
#
  sed s/\'kpp\'/\'kpp_mom4p1\'/g -i input.nml 
  sed "s/ocean_vert_kpp_nml/ocean_vert_kpp_mom4p1_nml/g" -i input.nml 
#
#Do not specify the layout so that the model could run on any npes
#  sed "/.*layout.*/d" -i input.nml
