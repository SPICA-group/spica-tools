ENM
===============
Usage   
--------
.. parsed-literal::

    cg_spica ENM :strong:`CG_PDB` :strong:`CG_TOP` [:strong:`-maxr` :emphasis:`<float>`] [:strong:`-kENM` :emphasis:`<float>`] [:strong:`-pspica`] [:strong:`-v1`] [:strong:`-dssp` :emphasis:`<string>`] [:strong:`-aapdb` :emphasis:`<string>`]

Description
-----------
``ENM`` applies the elastic network model (ENM) to CG configuration and generates
topology files to prepare input files for running CG-MD with SPICA.
ENM requires values for cutoff length and force constants, 
specified with the ``-maxr`` and ``-kENM`` options, respectively. 
The default values are set to 9.0 A and 1.195 kcal/A\ :sup:`2`\/mol, and these values are used 
in the SPICA protein model (See the `SPICA protein ver1`_ or `SPICA protein ver2`_ paper).

The SPICA protein ver1 does not have backbone parameters dependent on protein secondary structure.
So, when using the ver1 parameter set, add an option ``-v1`` for not applying DSSP.

For SPICA protein ver2, incorporating secondary structure-dependent backbone parameters, 
the `DSSP`_ program and protein atomistic configuration (in PDB format)
are needed. ``-dssp`` and ``-aapdb`` options take a file path to DSSP and the protein PDB, respectively.
When DSSP can be used through a command on a console, the ``-dssp`` option can be skipped.


.. _SPICA protein ver1: https://pubs.acs.org/doi/10.1021/acs.jctc.1c01207
.. _SPICA protein ver2: https://pubs.acs.org/doi/10.1021/acs.jctc.3c01016
.. _DSSP: https://swift.cmbi.umcn.nl/gv/dssp

Example
-------
.. code-block:: bash

    cg_spica ENM 2mag.cg.pdb 2mag.cg.top.v1 -v1

.. code-block:: bash

    cg_spica ENM 2mag.cg.pdb 2mag.cg.top.v2 -aapdb 2mag.aa.pdb -dssp /usr/bin/dssp

.. code-block:: bash

    cg_spica ENM 2mag.cg.pdb 2mag.cg.top.pspica -pspica

.. image:: figures/ENM_example.jpg
    :scale: 30
    
protein CG configuration (CPK) applied ENM (red)

:download:`2mag.cg.pdb <data/2mag.cg.pdb>` 
:download:`2mag.aa.pdb <data/2mag.aa.pdb>`  

:download:`2mag.cg.top.v1 <data/2mag.cg.top.v1>`  
:download:`2mag.cg.top.v2 <data/2mag.cg.top.v2>`  
:download:`2mag.cg.top.pspica <data/2mag.cg.top.pspica>`  

Positional args
---------------

``CG_PDB`` [<.pdb>] 
    Input protein CG configuration
``CG_TOP`` [<.top>] 
    Output SPICA topology

Optional args
-------------

``-maxr`` <float> (9.0) [A]
    cutoff length for applying ENM
``-kENM`` <float> (1.195) [kcal/A\ :sup:`2`\]
    force constant of harmonic potential for ENM
``-pspica`` (off)
    Assign partial charge (0.5590) for pSPICA FF (default: 0.1118, for SPICA FF)
``-v1`` (off)
    Generate topology file for SPICA ver1, not considering backbone parameters dependent on protein secondary structure
``-dssp`` <string>
    File path to the DSSP program
``-aapdb`` <string>
    Atomistic protein configuration for applying DSSP

