setup_gmx
===============
Usage   
--------
.. parsed-literal::

    cg_spica setup_gmx <:strong:`topfile 1`> <:strong:`nmol 1`> [<:strong:`topfile 2`> <:strong:`nmol 2`> ... <:strong:`topfile n`> <:strong:`nmol n`>] <:strong:`paramfile`> [<:strong:`coordinate`>] [:strong:`-prot`]

Description
-----------
``setup_gmx`` generates simulation input files for CG-MD with SPICA using GROMACS.
This program requires the following files as inputs:

- ``topfile``: SPICA topology file, distributed at `SPICA web <https://www.spica-ff.org/download.html>`_
- ``paramfile``: SPICA parameter file, distributed at `SPICA web <https://www.spica-ff.org/download.html>`_

The ``topfile`` can be generated with the ``json2top`` or ``ENM`` (for proteins only) commands.
For more details, see :doc:`json2top <json2top>` and :doc:`ENM <ENM>`.
The ``nmol`` is the number of ``topfile`` molecules in the system.  
It is sensitive to the order of arguments, so the molecule topology and number 
should be given in the same order as given in the prepared configuration file. 

Unlike the ``setup_lmp`` command, this program does not normally require a pdb 
configuration file, because it generates input files needed only to execute
the `gmx grompp <https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html>`_ 
command in GROMACS. 

Use the ``-prot`` option is used to set up systems containing proteins.
This is because SPICA takes protein equilibrium angle values from the initial configuration.
In this case, a configuration file in PDB format must be specified as the last 
argument of the program. 

Running the program with the appropriate input files yields the following outputs:

- ``toppar/``

  - ``SPICA.itp``: SPICA force field parameter
  - ``{molecule}.itp``: Included topology files for each molecule.

- ``topol.top``: System topology file including ``SPICA.itp`` and ``{molecule}.itp``
- ``CGindex.ndx``: Atom index file of a target system
- ``out.psf``: Topology file of a target system in PSF format

.. topic:: NOTE

    The standard GROMACS package does NOT support the SPICA angle or nonbonded interactions. 
    To use SPICA with GROMACS, the package must be modified with the patch file for the angular
    interaction in `gromacs-spica <https://github.com/SPICA-group/gromacs-spica>`_ repository,
    and the tabulated potentials for the nonbonded interaction must be applied.

Example
-------
.. code-block:: bash

    cg_spica setup_gmx DOPC.top 128 WAT.top 1408 spica_db.prm

:download:`DOPC.top <data/DOPC.top>`  
:download:`WAT.top <data/WAT.top>`  
:download:`spica_db.prm <data/spica_db.prm>`  

:download:`topol.top <data/topol.top>`  
:download:`SPICA.itp <data/toppar/SPICA.itp>`  
:download:`DOPC.itp <data/toppar/DOPC.itp>`  
:download:`WAT.itp <data/toppar/WAT.itp>`  
:download:`CGindex.ndx <data/CGindex.ndx>`  
:download:`out.psf <data/out.psf>`  

.. code-block:: bash

    cg_spica setup_gmx 1d6x.cg.top 1 DOPC.top 128 WAT.top 2134 CLA.top 4 spica_db.prm prot_memb.cg.pdb -prot

:download:`1d6x.cg.top <data/setup_gmx_prot/1d6x.cg.top>`  
:download:`DOPC.top <data/DOPC.top>`  
:download:`WAT.top <data/WAT.top>`  
:download:`CLA.top <data/CLA.top>`  
:download:`spica_db.prm <data/spica_db.prm>`  
:download:`prot_memb.cg.pdb <data/setup_gmx_prot/prot_memb.cg.pdb>`  

:download:`topol.top <data/setup_gmx_prot/topol.top>`  
:download:`SPICA.itp <data/setup_gmx_prot/toppar/SPICA.itp>`  
:download:`1d6x.cg.itp <data/setup_gmx_prot/toppar/1d6x.cg.itp>`  
:download:`DOPC.itp <data/setup_gmx_prot/toppar/DOPC.itp>`  
:download:`WAT.itp <data/setup_gmx_prot/toppar/WAT.itp>`  
:download:`CLA.itp <data/setup_gmx_prot/toppar/CLA.itp>`  
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

