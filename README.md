# spica-tools-devel
SPICA tools  
## Tools
* C  
  * setup_lammps  
  * setup_gromacs  
* Python (python2)
  * common  
    * map_to_cg  
    * gen_top_elastic_network  
      * usage: gen_top_elastic_network.py [-h] [-maxr MAXR] [-kENM KENM] [-pspica] input output  
         * positional arguments:  
           * input    Specify input CG PDB file name.  
           * output    Specify output topology file name.  
         * optional arguments:  
           * -h, --help  show this help message and exit  
           * -maxr MAXR  Cutoff length of ENM (default: 9.0 A).  
           * -kENM KENM  Force constant for ENM (default: 1.195 kcal/A2).  
           * -pspica     Assign partial charge (0.5990) for pSPICA FF (default: 0.1118, for SPICA FF).  
      * example: gen_elastic_network.py protein.cg.pdb protein.cg.top
    
  * pspica
    * WAT2PWAT  
    * convert_PARM    
