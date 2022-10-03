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

The ``topfile`` can be generated with ``json2top``, or ``ENM`` (for proteins only) commands.
For more information on these, see :doc:`json2top <json2top>` and :doc:`ENM <ENM>`.
The ``nmol`` is the number of ``topfile`` molecules included in your system.  
It is sensitive to the order of arguments, so the molecule topology and number 
should be given in the same order as given in the prepared configuration file in the PDB file. 

Running the program with the appropriate input files yields the following outputs:

- ``DATA.FILE``: LAMMPS DATA file 
- ``PARM.FILE``: Interaction parameter file that should be read when executing LAMMPS CG-MD. 
- ``out.psf``: Topology file of a target system in PSF format

``DATA.FILE`` and ``PARM.FILE`` are required to run LAMMPS CG-MD with SPICA FF. 
Specify them in the LAMMPS simulation input file to start the simulation. 
Input examples are distributed in
`lipid membrane tutorial <https://www.spica-ff.org/tutorial_lipid3.html>`_ on SPICA web.
Although ``out.psf`` is not needed to run CG-MD, it will be useful to visualization and
analysis of the system because it has all the topology information on the bonds and angles 
of molecules included in the system.

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


