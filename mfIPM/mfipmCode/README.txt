mfipm: Matrix-free Interior Point Method for Compressed Sensing problems.

[1] Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich, 
    "Matrix-free Interior Point Method for Compressed Sensing Problems",
    ERGO Technical Report, 2012.
---------------------------------------------------------------------------


1. Introduction
===============

The present software package contains the Matrix-free IPM solver which can be 
used for the solution of one dimensional Compressed Sensing problems. 
It implements the l1-regularised program:

  min tau*||x||_1 + 0.5*||Ax - b||^2.

Examples are given on how to use the Matrix-free IPM.

NOTE: All the MATLAB .m files support help info, for the help command of MATLAB.


2. Installation and Setup
=========================

Start Matlab and make sure that the working directory is set to the
main directory of the present package.  At the MATLAB prompt, run

  >> setupmfipm

This script adds the appropriate directories in the MATLAB path.  
The script will try to permanently add these directories to your 
path (in pathdefs.m), but may fail if that file is read-only.  
In that case, please copy and paste to your startup.m file the 
'addpath' commands printed to the screen.

The startup.m file is not by default created in the startup folder of 
MATLAB. It has to be created by the user. For more information see:
http://www.mathworks.co.uk/help/matlab/ref/startup.html.
To find your current MATLAB startup folder type in MATLAB prompt:
  >> userpath
the returned value will be the directory of the startup folder.
For information about MATLAB startup folder see:
http://www.mathworks.co.uk/help/matlab/matlab_env/matlab-startup-folder.html.

To test if the installation and setup for the matrix-free ipm have been 
completed successfully please run in the MATLAB prompt:

  >> mfipm_demo


3. Quick Start
==============

The standard call of the Matrix-free IPM is:

  [x out] = mfipm(n, A, b, tau, opts) 

The first four arguments (n, A, b and tau) are required. The opts argument
is a structure and it is optional. Here is a short explanation of the available 
options for opts:

  opts.tol:        Optimality tolerance of the Matrix-free IPM.
                   Default value: 1.0e-8.
  opts.maxiters:   Maximum iterations of the Matrix-free IPM.
                   Default value: 100.
  opts.tolpcg:     Optimality tolerance for the PCG method.
                   Default value: 1.0e-4.
  opts.mxiterpcg:  Maximum iterations of the PCG method.
                   Default value: 200.
  opts.maxstepsz:  Maximum stepsize for primal and dual directions.
  opts.cnt1:       Centering parameter when predictor directions are 
                   performed. Set to 0 for pure predictor directions.
                   Default value 0.1;
  opts.cnt2:       Centering parameter when predictor directions are 
                   performed and one of the primal and dual stepsizes of 
                   the previous iteration were less than opts.ap for last
                   iteration. 
                   Set to 0 for pure predictor directions.
                   Default value 0.5;
  opts.cnt3:       Centering parameter when corrector directions are 
                   performed. Set to 1 for perfect centering to the
                   central path.
                   Default value 0.8;
  opts.ap:         If primal or dual stepsize <= ap then use cnt2 as as
                   centering parameter, otherwise use cnt1.
                   Default value 0.5;
  opts.ac:         If primal or dual stepsize after calculating the
                   predictor direction are <= ac then perform a corrector
                   update.
                   Default value 0.1;
  opts.linsol:     0 to use the Conjugate Gradients method as a linear 
                   system solver, 1 to use a direct solver.
                   Default value: 0.
  opts.verbose:    0 for no printing info, 1 for basic printing info
                   and 2 for detailed printing info.
                   Default value: 1.

The list with the output parameters follows.

  x:               The reconstructed sparse representation.
  out.totcputime:  Total CPU time required by the Matrix-free IPM.
  out.iters:       Number of Matrix-free IPM iterations.
  out.cpusteps:    CPU time used for the calculation of the stepsizes.
  out.totpcgiters: Total number of PCG iterations.
  out.nMat:        Number of matrix-vector products.
  out.predictors:  Number of predictor directions.
  out.correctors:  Number of corrector directions.
  out.linsystems:  Total number of linear systems solved. It is equal
                   to out.predictors + out.correctors.
  out.status:      1 if the problem was solved successfully, 0 if 
                   maximum number of iterations reached.
  out.message:     Termination status of the Matrix-free IPM
                   solver. One of the following two messages is returned:
                   i)  'Converged'
                   ii) 'Maximum iterations reached'

The user should have a look at the "mfipm_demo.m". In this file two 
simple cases of how to use the matrix-free ipm are demonstrated.


4. External packages
===========================================

Directory "operators"
   * Subdirectory "@A_operator":
     A wrapper for matrix-operators in order to be used as implicit matrices.
     Copied from the FPC_AS_v1.21 package.
     url: http://www.caam.rice.edu/~optimization/L1/FPC_AS/
   * Subdirectory "PartialDCT":
        > "pdct.m": partial 1-d discrete cosine transform operator.
     	   Copied from the FPC_AS_v1.21 package.
     	   url: http://www.caam.rice.edu/~optimization/L1/FPC_AS/


5. License
==========

 MFIPM or Matrix-free Interior Point Method
 Copyright (C) 2013, Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.


6. The Authors
==============

If you have any bug reports or comments, please feel free to email one 
of the authors:

  Kimon Fountoulakis <K.Fountoulakis@sms.ed.ac.uk>
  Jacek Gondzio      <J.Gondzio@ed.ac.uk>
  Pavel Zhlobich     <P.Zhlobich@ed.ac.uk>

Enjoy!
Kimon, Jacek and Pavel
17 March 2013
