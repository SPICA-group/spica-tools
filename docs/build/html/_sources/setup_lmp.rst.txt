setup_lmp
===============
Usage   
--------
.. parsed-literal::

    cg_spica setup_lmp <:strong:`topfile 1`> <:strong:`nmol 1`> [<:strong:`topfile 2`> <:strong:`nmol 2`> ... <:strong:`topfile n`> <:strong:`nmol n`>] 
                       <:strong:`paramfile`> <:strong:`coordfile`>

Description
-----------
``setup_lmp`` generates simulation input files for CG-MD with SPICA using LAMMPS.
This program requires the following files as inputs:

- ``topfile``: SPICA topology file, distributed at `SPICA web <https://www.spica-ff.org/download.html>`_
- ``paramfile``: SPICA parameter file, distributed at `SPICA web <https://www.spica-ff.org/download.html>`_
- ``coordfile``: CG configuration file in PDB format

``topfile`` can be generated with ``json2top``, or ``ENM`` (only for proteins) commands.
For details about them, see :doc:`json2top <json2top>` and :doc:`ENM <ENM>`.
``nmol`` is the number of ``topfile`` molecules included in your system.  
**Note** that one should be careful with the order of arguments. Especially, molecule 
topology and number should be the same as in the PDB file. 

Executing the program with proper input files, you can obtain the following outputs:

- ``DATA.FILE``: LAMMPS DATA file 
- ``PARM.FILE``: Interaction parameter file, which should be read when performing LAMMPS CG-MD. 
- ``out.psf``: Topology file of a target system in PSF format

``DATA.FILE`` and ``PARM.FILE`` will be mandatory to run LAMMPS CG-MD with SPICA FF. They are specified in
a LAMMPS simulation input file to launch simulations. Example inputs are distributed at 
`lipid membrane tutorial <https://www.spica-ff.org/tutorial_lipid3.html>`_ in SPICA web.
``out.psf`` is not needed in terms of running CG-MD, but it would be useful to visualize or analyze your 
system because it has all topology information on bonds and angles of molecules included in the system.

Example
-------
``cg_spica setup_lmp DOPC.top 128 WAT.top 1408 spica_db.prm dopc.cg.pdb`` 

:download:`DOPC.top <data/DOPC.top>`  
:download:`WAT.top <data/WAT.top>`  
:download:`spica_db.prm <data/spica_db.prm>`  
:download:`dopc.cg.pdb <data/dopc.cg.pdb>`  

:download:`DATA.FILE <data/DATA.FILE>`  
:download:`PARM.FILE <data/PARM.FILE>`  
:download:`out.psf <data/out.psf>`  

Positional args
---------------

``topfile`` <topology>
    topology file of a molecule
``nmol`` <int>
    number of molecules
``paramfile`` <parameter>
    SPICA force field parameter file
``coordfile`` <.pdb>
    CG configuration file


