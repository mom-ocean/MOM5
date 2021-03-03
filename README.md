# The Modular Ocean Model + TIDES

MOM is a numerical ocean model based on the hydrostatic primitive equations. Development of the model is managed through the Model Development Lab at
[http://www.mom-ocean.org](http://www.mom-ocean.org)

Explicit tidal forcing in OGCM (ocean general circulation models) using  tides in MOM5 plus accurate phase information is introduced. 


Some tips to use tides in MOM5

* must use  vertical_coordinate z* to enable tides in the surface
* turn off Bryan Lewis vertical diffusity in ocean_vert_mix 


## 8 dominant constituents of semidiurnal and diurnal tides plus phase information
[ocean_barotropic.F90](https://github.com/mabelcalim/mom/blob/master/src/mom5/ocean_core/ocean_barotropic.F90) with phase

[input.nml](https://github.com/mabelcalim/mom/blob/master/nml/CGCM_T8/input.nml)
[data_table](https://github.com/mabelcalim/mom/blob/master/nml/CGCM_T8/data_table)
[diag_table](https://github.com/mabelcalim/mom/blob/master/nml/CGCM_T8/diag_table)
[field_table](https://github.com/mabelcalim/mom/blob/master/nml/CGCM_T8/field_table)


## 11 tidal constituents = 8 dominant + 3 long term (Ssa, Mf and Mm)
[ocean_barotropic_longterm.F90](https://github.com/mabelcalim/mom/blob/master/src/mom5/ocean_core/ocean_barotropic_longterm.F90)



## Using tide from a file
Regriding tideamp.nc file using TPX7.2 OSU Tidal Data Inversion form Oregon as a example step by step in this [link](https://nbviewer.jupyter.org/github/mabelcalim/notebooks/blob/master/tideamp.nc.ipynb?create=1)!



# More info 
My thesis [PT] : ["The global barotropic tidal-driven impact on climate timescales"](http://mtc-m21b.sid.inpe.br/col/sid.inpe.br/mtc-m21b/2017/06.05.17.39/doc/publicacao-1.pdf) 





