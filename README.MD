<div align="center">
<img src="https://raw.githubusercontent.com/tpms-lattice/ASLI/master/docs/manual/figures/ASLI.png" height="142">
</div>

<p> <h3 align="center"> ASLI (A Simple Lattice Infiller) is an open source tool for creating lattice infills with varying unit cell type, size and/or feature</h3> &nbsp; </p>

<div align="center" markdown="1">

[![Website][website-image]][website] [![Release][release-image]][releases] [![License][license-image]][license]

[release-image]: https://img.shields.io/github/v/release/tpms-lattice/ASLI?color=blue&label=release&style=flat-square
[releases]: https://github.com/tpms-lattice/ASLI/releases

[license-image]: https://img.shields.io/badge/license-AGPL-blue.svg?style=flat-square
[license]: https://github.com/tpms-lattice/ASLI/blob/master/LICENSE

[website-image]: https://img.shields.io/static/v1?label=ASLI&message=website&color=blue&style=flat-square
[website]: http://www.biomech.ulg.ac.be/ASLI
</div>

<div align="center">
  <img src="https://img.shields.io/static/v1?label=Ubuntu&logo=Ubuntu&logoColor=white&message=support&color=success&style=flat-square" alt="Ubuntu support">
  <img src="https://img.shields.io/static/v1?label=Windows&logo=windows&logoColor=white&message=support&color=success&style=flat-square" alt="Windows support">
</div>

---

## Getting started with ASLI ##
The easiest way to get started with ASLI is to use the pre-build binaries found in [releases](https://github.com/tpms-lattice/ASLI/releases). For more advanced advanced users, it may be more interesting to build ASLI from scratch. The prerequisites and procedure to build ASLI are detailed below.

### Prerequisites ###
ASLI is written in C++ and makes use of [CMake](http://www.cmake.org) and [GNU Make](https://www.gnu.org/software/make/) to automate the building process. ASLI also has the following dependencies:
  * [AdaptTools](https://github.com/ISCDtoolbox/AdaptTools)
  * [ALGLIB](https://www.alglib.net)
  * [CGAL](https://www.cgal.org)
    * [GMP](https://gmplib.org/)
    * [MPFR](https://www.mpfr.org/)
    * [Boost](https://www.boost.org/)
    * [TBB](https://intel.com/oneTBB) (optional, necessary for parallelization)
  * [Eigen](http://eigen.tuxfamily.org)
  * [MMG](https://www.mmgtools.org)
  * [TETGEN](http://tetgen.org)
  * [Yaml-cpp](https://github.com/jbeder/yaml-cpp)
  * [QASLI](https://github.com/tpms-lattice/QASLI) (optional, the GUI of ASLI)
    * [QT](https://www.qt.io/) (optional, required if compiling the GUI)

Most of these dependencies are included with ASLI so that the user does not need to worry about them. Not included with ASLI are GMP, MPFR, Boost and TBB (optional). QT, required to compile the GUI of ASLI is not included either. To be able to compile ASLI, depending on the compilation settings, users will need some or all of these libraries available on their system.

### Install dependencies ###
#### Linux ####
If you are missing some required or optional dependencies they can be installed on Debian and Ubuntu with a couple of commands:

Install build tools:  
```Shell
sudo apt-get install git cmake
```

Install compilers:  
```Shell
sudo apt-get install g++
```

Install libraries (GMP, MPFR, Boost):  
```Shell
sudo apt-get install libgmp-dev libmpfr-dev libboost-all-dev
```

Install libraries (TBB):  
```Shell
sudo apt-get install libtbb-dev
```

Install libraries (QT):  
```Shell
sudo apt-get install qt3d5-dev
```

#### Windows (MSYS2) ####
In windows compilation is supported with MSYS2. If you are missing some required or optional dependencies they can be installed in MSYS2 with a couple of commands:

Install build tools:  
```Shell
pacman -S git cmake mingw-w64-x86_64-cmake
```

Install compilers, and libraries (GMP, MPFR, Boost):  
```Shell
pacman -S --needed base-devel mingw-w64-x86_64-toolchain
pacman -S msys2-runtime-devel
```

Install libraries (TBB):  
```Shell
pacman -S mingw-w64-x86_64-intel-tbb
```

Install libraries (QT):  
```Shell
pacman -S mingw-w64-x86_64-qt5
```

Add the location of tbb.dll and tbbmalloc.tbb to the windows environment variables.

### Downloading and compiling ASLI ###
#### Linux ####
Assuming [CMake](http://www.cmake.org) as well as [GNU Make](https://www.gnu.org/software/make/) are installed on the system and all required dependencies are available the main steps to build are:

1. Retrieve ASLI from the repository:  
    ```Shell
    git clone https://github.com/tpms-lattice/ASLI.git
    ```

2. Compile:  
    ```Shell
    cd ASLI
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DMARCH_NATIVE=ON ..
    make
    ```

To clean up, simply delete the `build` and `bin` directory found in the source directory.

#### Windows ####
Assuming [MSYS2](https://www.msys2.org/) is installed on the system and all required dependencies are available the main steps to build are:

1. Retrieve ASLI from the repository:  
    ```Shell
    git clone https://github.com/tpms-lattice/ASLI.git
    ```

2. Compile:  
    ```Shell
    cd ASLI
    mkdir build
    cd build
    cmake -G "MSYS Makefiles" -DCMAKE_BUILD_TYPE=Release -DMARCH_NATIVE=ON ..
    make
    ```

To clean up, simply delete the `build` and `bin` directory found in the source directory.

### Parallel mode ###
To compile ASLI in parallel mode include `-DCGAL_ACTIVATE_CONCURRENT_MESH_3=ON` when calling cmake. Parallel mode is currently limited to the CGAL workflow.

### Graphical user interface ###
To compile ASLI together with its GUI include `-DASLI_GUI=ON` when calling cmake.

## Running ASLI ##
To run ASLI from the command line you only need to type ```./ASLI config.yml``` if using Linux or ```ASLI.exe config.yml``` if using windows. The config.yml file specified points to the configuration file of ASLI.

## Issues and questions ##
If you find a bug please post an [issue](https://github.com/tpms-lattice/ASLI/issues). If you have questions regarding ASLI you can post them on [stackoverflow](https://stackoverflow.com/questions/tagged/ASLI) with the tag [`ASLI`](https://stackoverflow.com/questions/tagged/ASLI).

## Contact ##
Contact: asli (at) kuleuven . be

## Reference ##
Perez-Boerema, F., M. Barzegari and L. Geris (2022). A flexible and easy-to-use open-source tool for building functionally graded scaffolds, *Virtual and Physical Prototyping*, 17:3, 682-699, DOI: [10.1080/17452759.2022.2048956](https://doi.org/10.1080/17452759.2022.2048956).

## License ##
ASLI is licensed under the terms of the [GNU Affero General Public License ](https://www.gnu.org/licenses).

Copyright © KU Leuven 2019-2022