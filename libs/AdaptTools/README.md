# AdaptTools [![Build Status](https://travis-ci.org/ISCDtoolbox/AdaptTools.svg?branch=master)](https://travis-ci.org/ISCDtoolbox/AdaptTools)

This directory contains a set of tools for mesh adaptation in two and three dimensions. They can be used and combined in a Linux/Unix/Posix-compliant shell script wrapper or called as functions from C/C++/F programs.

#### Mesh adaptation
Mesh adaptation strategy has proved to be very efficient in numerical solutions of PDEs. The ability to adapt the elements geometric properties (size, shape and orientation) according to specific quantities related to the data greatly improves the accuracy of the solution and the convergence rate of numerical schemes.

The approach pionneered in our group (F. Hecht et al.) relies on the construction of nearly uniform simplicial meshes in a metric space. The underlying metric is related to the (approximate) Hessian of the solution of the problem.
![aniso](https://cloud.githubusercontent.com/assets/11386916/13727601/ecdd2862-e8f8-11e5-97ed-1e436c98a851.png)

#### Installation
Download the zip archive of AdaptTools or clone this repository:

   ` git clone https://github.com/ISCDtoolbox/AdaptTools.git `

navigate to the downloaded directory: 

   ` cd AdaptTools `

create a build directory and compile with cmake

   ```
   mkdir build
   cd build
   cmake ..
   make
   make install
   ```

#### Authors & contributors
* these tools have been developed and contributed over the years by several contributors. Main developers are Charles Dapogny (Université J. Fourier), Cécile Dobrzynski (Université de Bordeaux, INRIA) and Pascal Frey (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
AdaptTools is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).

