.. spica-tools documentation master file, created by
   sphinx-quickstart on Tue Sep 13 18:18:02 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to spica-tools' documentation!
=======================================
Spica-tools is a suite of programs that support the setup and parameterization of
coarse-grained (CG) molecular dynamics (MD) simulations using SPICA force field 
(Surface Property fItting Coarse grAined model, see `SPICA-web <https://www.spica-ff.org/>`_ for details). 
The tools can be used with the ``cg_spica`` command integrated python codes.
Spica-tools github repository is found at 
`spica-tools <https://github.com/SPICA-group/spica-tools>`_.


Installation           
=======================================
The command to invoke the tools in the src/ directory is cg_spica and requires 
the Python 3 environment. 
You will need a bash environment to execute the tools using the command in a console. 
You can download this repository, including the setup file to use the command,
with the following command:

``git clone git@github.com:SPICA-group/spica-tools.git``

Then, change the current directory to the top directory of the repository:

``cd spica-tools/``  

There is a bash script set.sh, to add the path to use the command, 
run it with the following command:

``source set.sh``  

The path to cg_spica has been added to the PATH environmental variable. 
Adding the above line to a configuration file such as files *.bash_profile* may
be useful if you need to use the spica-tools frequently. 
To verify that the commands to integrate python codes in the src/ directory 
is available in the console, run the following commands:

``cg_spica -h``  

``cg_spica map2cg -h``  

You can see the available options for the commands, if the command path 
has been correctly added to the PATH. Cg_spica requires several python libraries to 
execute some codes, including arithmetic operations, and to read trajectory files 
obtained from MD simulations. These libraries are listed in requirements.txt in 
the top directory and can be installed using pip with the following command:

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
