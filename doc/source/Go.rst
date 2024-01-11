Go
===============
Usage   
--------
.. parsed-literal::

    cg_spica Go :strong:`CG_PDB` :strong:`CG_TOP` [:strong:`-maxr` :emphasis:`<float>`] [:strong:`-eps` :emphasis:`<float>`] [:strong:`-pspica`] [:strong:`-v1`] [:strong:`-dssp` :emphasis:`<string>`] [:strong:`-aapdb` :emphasis:`<string>`]

Description
-----------
``Go`` applies the Go-like model to protein CG configuration and generates
topology files to prepare input files for running CG-MD with SPICA.
Go requires values for cutoff length and interaction energy, 
specified with the ``-maxr`` and ``-eps`` options, respectively. 
The default values are set to 9.0 A and 1.5 kcal/mol. 
The other options such as ``-dssp`` and ``-aapdb`` are the same as the :doc:`ENM <ENM>` command.


Please note that the SPICA protein with the Go-like model has not been published yet.


Example
-------
.. code-block:: bash

    cg_spica Go 2mag.cg.pdb 2mag.cg.top.v1.Go -v1

.. code-block:: bash

    cg_spica Go 2mag.cg.pdb 2mag.cg.top.v2.Go -aapdb 2mag.aa.pdb
    
:download:`2mag.cg.pdb <data/2mag.cg.pdb>` 
:download:`2mag.aa.pdb <data/2mag.aa.pdb>` 

:download:`2mag.cg.top.v1.Go <data/2mag.cg.top.v1.Go>`  
:download:`2mag.cg.top.v2.Go <data/2mag.cg.top.v2.Go>`  


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
``-eps`` <float> (1.5) [kcal/mol]
    interaction energy of applied Go-like model
``-pspica`` (off)
    Assign partial charge (0.5590) for pSPICA FF (default: 0.1118, for SPICA FF)
``-v1`` (off)
    Generate topology file for SPICA ver1, not considering backbone parameters dependent on protein secondary structure
``-dssp`` <string>
    File path to the DSSP program
``-aapdb`` <string>
    Atomistic protein configuration for applying DSSP

