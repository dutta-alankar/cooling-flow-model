# cooling-flow-model-
This work is on modelling of local cooling flows in astrophysical plasma of the Circumgalactic medium and comparison with Illustris TNG50 Cosmological simulation.

```
   ___   ___     ___   __    __ __  __   ___      ____ __      ___   __    __
  //    // \\   // \\  ||    || ||\ ||  // \\    ||    ||     // \\  ||    ||
 ((    ((   )) ((   )) ||    || ||\\|| (( ___    ||==  ||    ((   )) \\ /\ //
  \\__  \\_//   \\_//  ||__| || || \||  \\_||    ||    ||__|  \\_//   \V/\V/ 
                                                                       
```

This repo contains the codes and the data used in the paper "Cooling flows around cold clouds in the circumgalactic medium:theory & comparison with TNG50"
The larger sized data used in the code and cannot be uploaded to Github due to size restrictions can be found here: https://drive.google.com/file/d/1a2xM4T_U-emlXu8vRK23LeGUMgklnoFw/view?usp=sharing
The ```cooling-flow-model-largedata.zip``` file contains the following files:
```
data-gasandcloud.h5
segmentation_67_h8_Mg-II-numdens_gt_1e-08.hdf5
tag-neighbours.npy
```
These files need to be placed at proper location for the python scripts to find them. You may need to tweak a little bit for the code to find it.

The Illustris TNG50 data that has been used is publicly hosted at: https://www.tng-project.org/

Description of the files:

- ```2-fluid-Riemann.nb```: A Mathematica notebook solving the Riemann problem of the one dimensional two-fluid model.
- ```bernoulli-model(rTBeMdot).txt```: Bernouli number for the solution of the cooling flow model.
- ```cooling-flow-MHD-model.py```: Generates cooling flow model and compares the model with stacked clouds from Nelson et. al 2020.
- ```cooltable.dat```: CLoudy generated cooling function used in the cooling flow modelling.
- ```discriminant_contour_MHD.py```: Contour of the discriminant of the two-fuid equation for the transonic solution to exist.
- ```discriminant_transonic_cf.py```: Pure hydro discriminant for the transonic solution to exist.
- ```pres_dr_CF_MHD-vary_for_figs1and3.py```: Creates just the de-dimensionalized solutions of fluid field in cooling flow model. Set beta very high to get the pure hydro limit. This can be varied to get the figures 1 and 3 in the paper.
- ```dylan-nelson/```: Contains cloud stacking data from Nelson et. al 2020.
- ```illustris-analysis/cloud-read.py```: Creates the essential data file ```data-gasandcloud.h5``` from the Illustris TNG501 snapshot. This creates a smaller file from which most analysis can be done without depending on the much larger snapshot. Requires Illustris data to run.
- ```illustris-analysis/diff-emm-comparepdfs.png```: Generates the normalized mass and volume PDFs and compares it with normalized differential emission from the entire halo.
- ```illustris-analysis/diff-emm-create_data.py```: MPI aware code which calculates the differential emission around a given region from the center of each cloud. Usage: ```mpiexec -n <# processors> python -m mpi4py diff-emm-create_data.py <cutoff in physical kpc>```. If ```<cutoff in physical kpc>``` is not provided, it uses 3 times the respective cloud radius as a default cutoff. The data is stored in a ```.npy``` file whose name is self-explanatory.
- ```illustris-analysis/diff-emm-plot_data.py```: PLots the differential emission data along with other useful stuff.
- ```illustris-analysis/diff-emm-neightag-nocloud.py```: MPI aware code that calculates the differential emission from regions of the halo which are not part of any cloud given some cutoff. Usage is similar to ```illustris-analysis/diff-emm-create_data.py```. Cells are also tagged and the data stored in a ```.npy``` file. Tag for cells are zero if it is part of no cloud and greater than zero depending upon how many overlapping clouds a cell in the halo is a part of.
- ```illustris-analysis/gamma_m-regress.py```: This illustrates the magnetic pressure and density relation for all the halo gas and does fitting of a power law to it 
for getting the polytropic index value for the magnetic fluid.
- ```illustris-analysis/gamma_m-regress-subsample.py```: Similar in functionality as ```illustris-analysis/gamma_m-regress.py``` but downsamples the showed data points for faster plot generation. Some fixed number of points are uniform randomly chosen in equally intervaled bins in the scatter plot and the color of each of the data points are smoothed by a Lowess filter.
- ```illustris-analysis/mdot-dn.py```: Calculates the mass inflow rate from the stacked cloud profiles used in Nelson et. al 2020.
- ```illustris-analysis/oneoff-cloud-DEM.py```: Creates the differential emission from the region chosen around the center of one user chosen cloud.
- ```illustris-analysis/power-spec.py```: Generates the gas velocity power spectrum of gas from the regions around the cloud. Useful for any future turbulence studies.
- ```illustris-analysis/proj-prof.py```: Generates projected histograms of the halo for a quick visual inspection. By all means use the plotting utilities from Illustris web interface and avoid this.
