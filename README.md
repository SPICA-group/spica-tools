# spica-tools-devel
This repository includes tools to set up coarse-grained (CG) molecular dynamics simulation 
with the SPICA force field and support SPICA CG parametrization.
The tools are basically written in Python and can be called through a command explaned below.
These tools were developed by members of SPICA-group. Details about the SPICA force field and 
the members can be found at (https://www.spica-ff.org).  

## Installation  
The command to call the tools in src directory is `cg_spica` and requires bash and python 3 enviroment.
You can dowload this repository including a setup file to use the command by the following command:  

    git clone git@github.com:SPICA-group/spica-tools-devel.git  
    
Then change the current directory to the top directory of the repository:  

    cd spica-tools-devel  
    
You can find a bash script `set.sh`, to add a path for using the command, and excute it by the following command:  

    source set.sh  

You can see 

## Tools
* Python 
  * json_to_top
  * map_to_cg  
  * map_traj
  * gen_top_ENM
  * WAT2PWAT  

* C  
  * setup_lammps  
  * setup_gromacs  
