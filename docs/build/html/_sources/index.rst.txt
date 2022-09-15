.. spica-tools documentation master file, created by
   sphinx-quickstart on Tue Sep 13 18:18:02 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to spica-tools' documentation!
=======================================
Spica-tools are to set up coarse-grained (CG) molecular dynamics (MD) simulation with 
SPICA (Surface Property fItting Coarse grAined model, see `SPICA-web <https://www.spica-ff.org/>`_
for details) and support the SPICA parametrization. 
The tools can be used with the ``cg_spica`` command integrating python codes.
The github repository of spica-tools is found at 
`spica-tools-devel <https://github.com/SPICA-group/spica-tools-devel>`_.


Installation           
=======================================
The command to call the tools in src/ directory is cg_spica and requires 
Python 3 enviroment. 
Bash enviroment is also needed to excute the tools with commands on consoles. 
You can dowload this repository including a setup file to use the command 
by the following command:

``git clone git@github.com:SPICA-group/spica-tools-devel.git``

Then change the current directory to the top directory of the repository:

``cd spica-tools-devel/``  

You can find a bash script set.sh, to add a path for using the command, 
and excute it by the following command:

``source set.sh``  

The path to cg_spica has been added to the environmental variable PATH. 
Adding the above line in your setup files like *.bash_profile* would 
be convinient if you often need to use the spica-tools. 
To make sure that the command integrating python codes in src/ directory 
is available on consoles, run the following commands:

``cg_spica -h``  

``cg_spica map2cg -h``  

You can see the available option of the commands by this if the command path 
is properly added to your PATH. Cg_spica requires several python libraries to 
excute some codes that include arithmetic operation and load trajectory files 
obtained from MD simulation. The libraries are listed in requirements.txt in 
the top directory and can be installed with pip by the following command:

``pip install -r requirements.txt``


Available command-line tools
=======================================
.. toctree::
   :maxdepth: 1

   json2top
   map2cg
   maptraj
   ENM
   wat2polar
   setup_lmp
   setup_gmx

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
