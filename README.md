# BRIC TUS Simulation Tools

MATLAB-based functions for running TUS acoustic simulations using T1-weighted MR and pseudo-CT images acquired at the Brain Research & Imaging Centre, University of Plymouth.

Currently, the simulation functions work for the NeuroFUS PRO CTX-500 and CTX-250 4-element transducers (https://brainbox-neuro.com/products/neurofus) only. The transcranial simulation function only supports input images with 1 mm isotropic voxels and the simulated pressure field is based on a free-field simulation at 20 W/cm<sup>2</sup> using the CTX-500 transducer. 


## Platform

Tested on macOS Catalina (10.15.7) and Monterey (12.5.1).


## Dependencies

* [k-Wave] (http://www.k-wave.org)
* [kArray tools] (http://www.k-wave.org/downloads/kWaveArray_alpha_0.3.zip)


## Instructions

Clone or download functions and utils and add them and the k-Wave and kArray tools to your MATLAB path.

### Free-field acoustic simulations
To run free-field acoustic simulations, run the following in MATLAB:
```
transducer = 'CTX500';
pressure = 51590;
phase = [0,319,278,237];
tussim_water_3D(transducer, pressure, phase)
```
**Input parameters:**
* `transducer`: CTX transducer type. This will determine the transducer central frequency and dimensions. Options are 'CTX500' or 'CTX250'.
* `pressure`: Source pressure applied by transducer in Pa.
* `phase`: 4-element array of phases of each transducer element in degrees for the focal depth required.

This function will run the free-field (i.e. through water) acoustic simulation for the input parameters you provide. The function will produce plots of the simulated acoustic pressure field and display the maximum pressure [MPa], mechanical index (MI), distance from the rear surface of the transducer [mm], I<sub>SPPA</sub> [W/cm<sup>2</sup>], and I<sub>SPTA</sub> [mW/cm<sup>2</sup>]. 

You will need to know the pressure and phases of each element of the transducer. The phase of each element can be obtained from the NeuroFUS PRO Transducer Power Output (TPO) unit. The pressure will determine the free-field I<sub>SPPA</sub> for the simulation. Values for a focal depth of 62 mm are given in the example above.

### Transcranial acoustic simulations
To run transcranial acoustic simulations on the example dataset provided, download the example dataset to the desired folder on your computer, then run the following in MATLAB, replacing `filepath` with the path to the folder where you saved the example data:
```
subj_id = 'sub-test01';
t1_filename = fullfile(filepath, 'sub-test01_t1w.nii');
ct_filename = fullfile(filepath, 'sub-test01_pct.nii');
output_dir = filepath;
focus_coords = [99, 161, 202];
bowl_coords = [90, 193, 262];
focus_depth = 60;
transducer = 'CTX500';
tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth, transducer)
```
**Parameters:**
* `subj_id`: ID of the subject you are running the simulation for.
* `t1_filename`: Full file path to the T1-weighted MR image.
* `ct_filename`: Full file path to the CT (or pseudo-CT) image.
* `output_dir`: Full path to the directory where you want your simulation output to be saved.
* `focus_coords`: 3-element array of voxel coordinates of the desired TUS focus. Add 1 if reading these off a viewer that uses zero-indexing (MATLAB indexing starts from 1).
* `bowl_coords`: 3-element array of voxel coordinates of the centre of the transducer **base**. Add 1 if reading these off a viewer that uses zero-indexing (MATLAB indexing starts from 1).
* `focus_depth`: Distance from transducer **face** to intended focus in mm, rounded to the nearest integer.
* `transducer`: CTX transducer type. This will determine the transducer central frequency and dimensions. Options are 'CTX500' or 'CTX250'.

This will load the example data and produce a 3D view of the skull model and transducer. The focus in this example is the dorsal anterior cingulate. 

Run the function again with the acoustic simulation flag turned on:
`tussim_skull_3D(subj_id, t1_filename, ct_filename, output_dir, focus_coords, bowl_coords, focus_depth, transducer, 'RunAcousticSim', true)`

This will produce plots of the simulated acoustic pressure field in each plane. The function will also display simulation parameters along with the coordinates of maximum pressure, distance from the rear surface of the transducer [mm], maximum pressure [MPa], mechanical index (MI), I<sub>SPPA</sub> [W/cm<sup>2</sup>], I<sub>SPTA</sub> [mW/cm<sup>2</sup>], pressure at intended focus [MPa], I<sub>SPPA</sub> at intended focus [W/cm<sup>2</sup>], and -6dB focal volume [mm<sup>3</sup>].
The plots, a .csv file and a .mat file containing the simulation output will be saved in the `output_dir`.

Please see function help for optional input arguments (e.g. to run thermal simulation).

To use your own data with these scripts, please see "Preparing input files" section below.


### Preparing input files
Currently, the transcranial simulation function only supports input images with 1 mm isotropic voxels and the output is based on a free-field simulation at 20 W/cm<sup>2</sup> using the CTX-500 transducer. 

You can provide your own source magnitude and transducer element phase information using the optional input.

You will need to supply a co-registered T1-weighted MR image and CT (or pseudo-CT) image for use in the simulations. The image matrix size can be of any size and will be automatically adjusted in the function to allow space for placing your transducer within the simulation grid.

It is preferable that you provide a T1-weighted MR image that has had noise outside the head masked out (e.g. the one used for pseudo-CT generation), so that the sensor is set to within this head mask. Otherwise, it will set the sensor to all voxels with intensity > 0 in the T1-weighted MR image. Alternatively, you can also use a brain extracted T1-weighted MRI, but this means you will not be able to simulate temperature rise at the skull interface when running the thermal simulation.

## Citing this work

If you use the functions in your own work, please consider citing the code (see "Cite this repository" in the About panel on the top right of this page), the paper below, and the relevant papers from the k-Wave toolbox.

The simulations as implemented in these functions are described in the following paper.

>    Siti N. Yaakub, Tristan White, Eric Kerfoot, Lennart Verhagen, Alexander Hammers, Elsa Fouragnan. Pseudo-CTs from T1-weighted MRI for planning of low-intensity transcranial focused ultrasound neuromodulation. (in preparation)

The simulation functions use the k-Wave Toolbox and kArray tools developed by Bradley Treeby and Ben Cox (University College London) and Jiri Jaros (Brno University of Technology). k-Wave and kArray tools can be downloaded for free from http://www.k-wave.org (see also Dependencies section) and described in the following papers. 

> B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave-fields," J. Biomed. Opt., vol. 15, no. 2, p. 021314, 2010.
>
> B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling nonlinear ultrasound propagation in heterogeneous media with power law absorption using a k-space pseudospectral method," J. Acoust. Soc. Am., vol. 131, no. 6, pp. 4324-4336, 2012.
>
> E. S. Wise, B. T. Cox, J. Jaros, B. E. Treeby, "Representing arbitrary acoustic source and sensor distributions in Fourier collocation methods," J. Acoust. Soc. Am., 146 (1), pp. 278-288, 2019.


## Questions/feedback
Feedback welcome at siti.yaakub@plymouth.ac.uk
