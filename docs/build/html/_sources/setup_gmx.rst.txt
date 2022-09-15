setup_gmx
===============
Usage   
--------
.. parsed-literal::

    cg_spica setup_gmx <:strong:`topfile 1`> <:strong:`nmol 1`> [<:strong:`topfile 2`> <:strong:`nmol 2`> ... <:strong:`topfile n`> <:strong:`nmol n`>] 
                       <:strong:`paramfile`> [<:strong:`coordfile`>] [:strong:`-prot`]

Description
-----------
``setup_gmx`` generates simulation input files for CG-MD with SPICA using GROMACS.
This program requires the following files as inputs:

- ``topfile``: SPICA topology file, distributed at `SPICA web <https://www.spica-ff.org/download.html>`_
- ``paramfile``: SPICA parameter file, distributed at `SPICA web <https://www.spica-ff.org/download.html>`_
- ``coordfile``: CG configuration file in PDB format

``topfile`` can be generated with ``json2top``, or ``ENM`` (only for proteins) commands.
For details about them, see :doc:`json2top <json2top>` and :doc:`ENM <ENM>`.
``nmol`` is the number of ``topfile`` molecules included in your system.  
**Note** that one should be careful with the order of arguments. Especially, molecule 
topology and number should be the same as in the system configuration file. 

Unlike the ``setup_lmp`` command, this program normally does not require pdb 
configuration file because it generates input files that are needed just to execute
the `gmx grompp <https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html>`_ 
command of GROMACS. 
The ``-prot`` option is used to setup systems including proteins because in SPICA
the proteins' equilibrium angle values are taken from the initial configuration.
In this case, the configuration file in PDB format must be specified as the last 
argument of the program. 

Executing the program with proper input files, you can obtain the following outputs:

- ``SPICA.itp``: SPICA force field parameter
- ``molecule.itp``: Molecule's information, such as atom types/names, residue names, and topology.
- ``topol.top``: System topology file including ``SPICA.itp`` and ``molecule.itp``
- ``CGindex.ndx``: Atom index file of a target system
- ``out.psf``: Topology file of a target system in PSF format

.. topic:: NOTE

    The standard GROMACS package does NOT support the SPICA angle and nonbonded interactions. 
    To use SPICA with GROMACS, one must modify the package with a patch file found in 
    `gromacs-SPICA <https://github.com/SPICA-group/gromacs-SPICA>`_ repository for 
    the angle interaction, and apply tabulated potentials for the nonbonded interaction.

Example
-------
``cg_spica setup_gmx DOPC.top 128 WAT.top 1408 spica_db.prm`` 

:download:`DOPC.top <data/DOPC.top>`  
:download:`WAT.top <data/WAT.top>`  
:download:`spica_db.prm <data/spica_db.prm>`  

:download:`SPICA.itp <data/SPICA.itp>`  
:download:`molecule.itp <data/molecule.itp>`  
:download:`topol.top <data/topol.top>`  
:download:`CGindex.ndx <data/CGindex.ndx>`  
:download:`out.psf <data/out.psf>`  

``cg_spica setup_gmx 1d6x.cg.top 1 DOPC.top 128 WAT.top 2134 CLA.top 4 spica_db.prm prot_memb.cg.pdb -prot``

:download:`1d6x.cg.top <data/setup_gmx_prot/1d6x.cg.top>`  
:download:`DOPC.top <data/DOPC.top>`  
:download:`WAT.top <data/WAT.top>`  
:download:`CLA.top <data/CLA.top>`  
:download:`spica_db.prm <data/spica_db.prm>`  
:download:`prot_memb.cg.pdb <data/setup_gmx_prot/prot_memb.cg.pdb>`  

:download:`SPICA.itp <data/setup_gmx_prot/SPICA.itp>`  
:download:`molecule.itp <data/setup_gmx_prot/molecule.itp>`  
:download:`topol.top <data/setup_gmx_prot/topol.top>`  
:download:`CGindex.ndx <data/setup_gmx_prot/CGindex.ndx>`  
:download:`out.psf <data/setup_gmx_prot/out.psf>`  

Positional args
---------------

``topfile`` <topology>
    topology file of a molecule
``nmol`` <int>
    number of molecules
``paramfile`` <parameter>
    SPICA force field parameter file

Optional args
-------------

``-prot`` (off)
    Read a pdb configuration file to extract reference angle for protein models
``coordfile`` <.pdb>
    CG configuration file, required when the ``-prot`` option is set


