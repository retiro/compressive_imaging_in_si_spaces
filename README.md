# Compressive Imaging in Shift-Invariant Spaces

### [Paper](https://arxiv.org/abs/2106.00404)

[Tin Vlašić](https://www.fer.unizg.hr/en/tin.vlasic), [Damir Seršić](https://www.fer.unizg.hr/en/damir.sersic) <br>
Faculty of Electrical Engineering and Computing, University of Zagreb, Croatia

This is the official implementation of the paper "Single-Pixel Compressive Imaging in Shift-Invariant Spaces via Exact Wavelet Frames."

## 1. Introduction
Intro + an image with the comparison of the conventional and our B-spline method

## 2. Setup
### mfIPM solver
This software package uses an external package of the Matrix-free IPM (mfIPM) solver for the compressive sensing minimization problem. The solver's source code and additional scripts can be found [here](https://www.maths.ed.ac.uk/ERGO/mfipmcs/). For convenience, we provide the solver's source code in the _mfIPM_ subdirectory.

To setup the mfIPM solver, start MATLAB and set the working directory to the main directory of the mfIPM package (...\compressive_imaging_in_si_spaces\mfIPM\mfipmCode). At the MATLAB prompt, run
```
>> setupmfipm
```
This script adds the appropriate directories in the MATLAB path. The script will try to permanently add these directories to your path (in pathdefs.m), but may fail if that file is read-only. In that case, please copy and paste to your startup.m file the 'addpath' commands printed to the screen. If that does not help, please read the README.txt file provided in the _mfipmCode_ subdirectory.

To test if the setup for the mfIPM has been completed successfully, please run in the MATLAB prompt:
```
>> mfipm_demo
```
### Fast WHT
If the mfIPM solver has been successfully installed, please set the working directory to the main directory of the software package (...\compressive_imaging_in_si_spaces). To use a fast realization of the Walsh-Hadamard transform (WHT), build a C program fastWht.c into a MEX file by running
```
>> mex fastWht.c
```
If you have problems building the fastWht.c program, please refer to the MathWorks [documentation](https://ch.mathworks.com/help/matlab/ref/mex.html) on buliding a MEX function.

You can test the fast WHT with a simple example code:
```
x = rand(1, 16);
y = fastWht(x);
xHat = fastWht(y); % the transform should reproduce x
```

### MATLAB toolboxes
The software package uses MATLAB signal processing and wavelet toolboxes. Instructions on installing a MATLAB toolbox can be found [here](https://www.mathworks.com/products/matlab/add-on-explorer.html).

## 3. Quick Start


## 4. Add-Ons
### Standard test images
The standard test images of 512 x 512 size that were used in the numerical experiments can be found in the _standard_test_images_ subdirectory.

### Reconstruction results and the λ parameters
In the _results_lambda_wht_tol-8-2_ subdirectory, there are reconstruction results in .mat files that correspond to various sparsity inducing biorthogonal wavelets. The reconstruction results are given in the terms of the PSNR and SSIM. These results correspond to the results for _bior2.2_ and _bior4.4_ wavelets given by the bar charts in Figures 4, 5, 6 and 7 of the paper. Additionally, the provided .mat files store optimal values of the regularization parameter λ, which were determined experimentally for every setting such that they lead to the best reconstruction results in terms of the PSNR.



## 4. License
```
Compressive Imaging in Shift-Invariant Spaces
Copyright (C) 2022,  Tin Vlasic and Damir Sersic

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
```

## 5. Citation
If you find our work useful in your research, please cite:
```
@inproceedings{vlasic2021cisispaces,
    author = {Vlasic, Tin
              and Sersic, Damir},
    title = {Single-Pixel Compressive Imaging in
             Shift-Invariant Spaces via Exact Wavelet Frames},
    booktitle = {arXiv},
    eprint={2106.00404},
    year={2021}
}
```

## 6. Contact
If you have any questions, please feel free to email the authors.
