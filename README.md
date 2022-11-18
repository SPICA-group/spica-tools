# spica-tools  
This repository provides tools to set up coarse-grained (CG) molecular dynamics (MD) simulation 
using the SPICA force field and tools to support SPICA CG parametrization.
The tools are basically written in Python and can be invoked with the commands described below.  
The tools were developed by members of the SPICA group. More information about the SPICA force field and 
the members can be found at (https://www.spica-ff.org).  

## Installation  
The command to invoke the tools in the src/ directory is `cg_spica` and requires the 
**Python 3** environment. You will need a bash environment to execute the tools using the commands in a console. 
You can download this repository, including the setup file to use the command, with the following command:  

    git clone git@github.com:SPICA-group/spica-tools.git  
    
Then, change the current directory to the top directory of the repository:  

    cd spica-tools/  
    
There is a bash script `set.sh` to add the path to use the command, and run it with the following 
command:  

    source set.sh  
    
The path to cg_spica has been added to the `PATH` environmental variable. Adding the above line to a 
configuration file such as `.bash_profile` may be useful if you need to use the spica-tools frequently. 
To verify that the commands to integrate python codes in the src/ directory is available in the console, 
run the following commands:

    cg_spica -h  
    cg_spica map2cg -h  
    
You can see the available options for the commands, if the command path is properly added to the PATH.  
Cg_spica requires several python libraries to execute some codes, including arithmetic operations, and 
to read trajectory files obtained from MD simulations. These libraries are listed in `requirements.txt` 
in the top directory and can be installed using `pip` with the following command:  

    pip install -r requirements.txt  
    
## Documentation  

The above spica-tools documentation can be found at https://spica-group.github.io/spica-tools. 
You can create the documentation on your machine with the following commands:  

    cd doc/  
    make html  


## Tools
* Python (in `src/`, can be used with `cg_spica`)
  * `json2top`  : create the top files needed for the setup_*** programs from a json file 
                  with SPICA mapping information  
  * `map2cg`    : map AA configuration to CG, available only in PDB format
  * `maptraj`   : map AA MD trajectories to CG, requires `MDAnalysis` module (https://www.mdanalysis.org)
  * `ENM`       : generate top files of protein with elastic networks
  * `Go`       : generate top files of protein with Go model
  * `wat2polar` : convert SPICA CG water to pSPICA polar CG water
  * `setup_lmp` : generate input files to run CG-MD with LAMMPS (https://www.lammps.org)
  * `setup_gmx` : generate input files to run CG-MD with GROMACS (https://www.gromacs.org).   
                  Apply the patch file distributed in the repository of SPICA-group 
                  (https://github.com/SPICA-group/gromacs-SPICA) to GROMACS software to make it useful 
                  for the MD with SPICA-FF

* C (in `c/`, cannot be used with cg_spica, must be compiled by the unix command `make` in each directory) 
  * setup_lammps  
  * setup_gromacs  
  
  The setup_lammps/gromacs codes have the same functionality as setup_lmp/gmx described above, 
  respectively. They may be useful for generating MD input files for very large systems.  

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
