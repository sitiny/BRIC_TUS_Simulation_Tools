# BRIC TUS Simulation Tools

MATLAB-based functions for running TUS acoustic simulations using T1-weighted MR and pseudo-CT images acquired at the Brain Research & Imaging Centre, University of Plymouth. To create pseudo-CT images, see https://github.com/sitiny/mr-to-pct.

Currently, the simulation functions work for the NeuroFUS PRO CTX-500 4-element transducer (https://brainbox-neuro.com/products/neurofus) only. The transcranial simulation function has been tested on input images with 1 mm isotropic voxels and the simulated pressure field is based on a free-field simulation at 20 W/cm<sup>2</sup> using the CTX-500-4 transducer. 


## Platform

Tested on Ubuntu 22.04.2 LTS (Jammy Jellyfish).


## Dependencies

* [k-Wave version 1.4] (http://www.k-wave.org)


## Instructions

Clone or download functions and utils and add them and the k-Wave v1.4 tools to your MATLAB path.

### Free-field acoustic simulations
To run free-field acoustic simulations, run the following in MATLAB:
```
pressure = 51590;
phase = [0,319,278,237];
tussim_water_3D(transducer, pressure, phase)
```
**Input parameters:**
* `pressure`: Source pressure applied by transducer in Pa.
* `phase`: 4-element array of phases of each transducer element in degrees for the focal depth required.

This function will run the free-field (i.e. through water) acoustic simulation for the input parameters you provide. The function will produce plots of the simulated acoustic pressure field and display the maximum pressure [MPa], mechanical index (MI), distance from the rear surface of the transducer [mm], I<sub>SPPA</sub> [W/cm<sup>2</sup>], and I<sub>SPTA</sub> [mW/cm<sup>2</sup>]. 

You will need the pressure and phases of each element of the transducer. The phase of each element can be obtained from the NeuroFUS PRO Transducer Power Output (TPO) unit. The pressure will determine the free-field I<sub>SPPA</sub> for the simulation. Values for a focal depth of 62 mm are given in the example above.

### Transcranial acoustic simulations
To run transcranial acoustic simulations on the example dataset provided, download the example data (T1-weighted MRI: https://osf.io/download/xhne5 and pseudo-CT: https://osf.io/download/fytwk) to the desired folder on your computer, then run the following in MATLAB, replacing `filepath` with the path to the folder where you saved the example data:
```
t1_filename = fullfile(filepath, 'sub-test01_t1w.nii');
ct_filename = fullfile(filepath, 'sub-test01_pct.nii');
output_dir = filepath;
focus_coords_in = [99, 161, 202];
bowl_coords_in = [90, 193, 262];
focus_depth = 60;
tussim_skull_3D(t1_filename, ct_filename, output_dir, focus_coords_in, bowl_coords_in, focus_depth)
```
**Parameters:**
* `t1_filename`: Full file path to the T1-weighted MR image.
* `ct_filename`: Full file path to the CT (or pseudo-CT) image.
* `output_dir`: Full path to the directory where you want your simulation output to be saved.
* `focus_coords`: 3-element array of voxel coordinates of the desired TUS focus. Add 1 if reading these off a viewer that uses zero-indexing (MATLAB indexing starts from 1).
* `bowl_coords`: 3-element array of voxel coordinates of the centre of the transducer **base**. Add 1 if reading these off a viewer that uses zero-indexing (MATLAB indexing starts from 1).
* `focus_depth`: Distance from transducer **face** to intended focus in mm, rounded to the nearest integer.

This will load the example data and produce a 3D view of the skull model and transducer. The focus in this example is the dorsal anterior cingulate. The visualisation is done at the default PPW = 3. We recommend doing transcranial simulations at PPW = 6.

Run the function again with the acoustic simulation flag turned on and PPW set to 6:
`tussim_skull_3D(t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth, transducer, 'PPW', 6, 'RunAcousticSim', true)`

This will produce plots of the simulated acoustic pressure field in each plane. The function will also display simulation parameters along with the coordinates of maximum pressure, distance from the rear surface of the transducer [mm], maximum pressure [MPa], mechanical index (MI), I<sub>SPPA</sub> [W/cm<sup>2</sup>], pressure at intended focus [MPa], I<sub>SPPA</sub> at intended focus [W/cm<sup>2</sup>], and -6dB focal volume [mm<sup>3</sup>]. The -6dB focal volume and the pressure field will be output as nifti files for overlay on the T1-weighted MRI.
The plots, output images and a .csv and .mat file containing the simulation output will be saved in the `output_dir`.

Please see function help for optional input arguments (e.g. to run thermal simulation).

To use your own data with these scripts, please see "Preparing input files" section below.


### Preparing input files
The output is based on a free-field simulation at 20 W/cm<sup>2</sup> using the CTX-500-4 transducer. 

You can provide your own source magnitude and transducer element phase information using the optional input. This currently only works for 4 element transducers, but the script can be edited to support 2-element transducers.

You will need to supply a co-registered T1-weighted MR image and CT (or pseudo-CT) image for use in the simulations. Preferably, your input images should have 1 mm isotropic voxels. The image matrix size can be of any size and will be automatically adjusted in the function to allow space for placing your transducer within the simulation grid.


## Citing this work

If you use the functions in your own work, please consider citing this repository [![DOI](https://zenodo.org/badge/528941098.svg)](https://zenodo.org/badge/latestdoi/528941098), the paper below, and the relevant papers from the k-Wave toolbox.

The simulations as implemented in these functions are described in the following paper.

>    Yaakub, S. N., White, T. A., Kerfoot, E., Verhagen, L., Hammers, A., & Fouragnan, E. F. (2023). Pseudo-CTs from T1-weighted MRI for planning of low-intensity transcranial focused ultrasound neuromodulation: an open-source tool. Brain Stimulation, 16(1), p75-78. https://doi.org/10.1016/j.brs.2023.01.838

The simulation functions use the k-Wave Toolbox and kArray tools developed by Bradley Treeby and Ben Cox (University College London) and Jiri Jaros (Brno University of Technology). k-Wave and kArray tools can be downloaded for free from http://www.k-wave.org (see also Dependencies section) and described in the following papers. 

> B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave-fields," J. Biomed. Opt., vol. 15, no. 2, p. 021314, 2010.
>
> B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling nonlinear ultrasound propagation in heterogeneous media with power law absorption using a k-space pseudospectral method," J. Acoust. Soc. Am., vol. 131, no. 6, pp. 4324-4336, 2012.
>
> E. S. Wise, B. T. Cox, J. Jaros, B. E. Treeby, "Representing arbitrary acoustic source and sensor distributions in Fourier collocation methods," J. Acoust. Soc. Am., 146 (1), pp. 278-288, 2019.


## Questions/feedback
Feedback welcome at siti.yaakub@plymouth.ac.uk
