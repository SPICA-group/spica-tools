# spica-tools-devel
This repository includes tools to set up coarse-grained (CG) molecular dynamics (MD) simulation 
with the SPICA force field and support SPICA CG parametrization.
The tools are basically written in Python and can be called through a command explaned below.  
The tools were developed by members of SPICA-group. Details about the SPICA force field and 
the members can be found at (https://www.spica-ff.org).  

## Installation  
The command to call the tools in src directory is `cg_spica` and requires bash and python 3 enviroment.
You can dowload this repository including a setup file to use the command by the following command:  

    git clone git@github.com:SPICA-group/spica-tools-devel.git  
    
Then change the current directory to the top directory of the repository:  

    cd spica-tools-devel  
    
You can find a bash script `set.sh`, to add a path for using the command, and excute it by the following 
command:  

    source set.sh  
    
The path to cg_spica has been added to the environmental variable `PATH`. To make sure that the command 
integrating
python codes in src/ directory is available on consoles, run the following commands:

    cg_spica -h  
    
    cg_spica map2cg -h  
    
You can see the available option of the command by this if the command path is properly added to your PATH.  
Cg_spica requires several python libraries to excute some codes that include arithmetic operation and load 
trajectory files obtained from MD simulation. The libraries are listed in `requirements.txt` in the top 
directory and can be installed with `pip` by the following command:  

    pip install -r requirements.txt  
    
The documentation of cg_spica is now in prepareation.


## Tools
* Python (can be used with `cg_spica`)
  * `json2top`  : make top files, required for setup_*** programs, from a json file having SPICA mapping 
                  information  
  * `map2cg`    : map AA configuration to CG, only PDB format is available
  * `maptraj`   : map AA MD trajectory to CG, `MDAnalysis` module (https://www.mdanalysis.org) is required
  * `ENM`       : generate top file of protein with elastic network
  * `wat2polar` : convert SPICA CG water to pSPICA polar CG water
  * `setup_lmp` : generate input files to run CG-MD with LAMMPS (https://www.lammps.org)
  * `setup_gmx` : generate input files to run CG-MD with GROMACS (https://www.gromacs.org) that is modified 
                  for SPICA-FF by applying a patch file distributed in another repository of SPICA-group 
                  (https://github.com/SPICA-group/gromacs-SPICA)

* C (cannot be used with cg_spica, needed to compile by the unix command `make` in each directory) 
  * setup_lammps.cpp 
  * setup_gromacs.cpp  
  
  The setup_lammps/gromacs codes have the same function as setup_lmp/gmx described above, 
  respectively. These may be useful to generate MD input files for extremely large systems.  

## License

Tools included in this repository are distributed under the MIT license.  

    Copyright (c) 2022 SPICA-group

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
