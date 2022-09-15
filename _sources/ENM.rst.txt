ENM
===============
Usage   
--------
.. parsed-literal::

    cg_spica ENM :strong:`CG_PDB` :strong:`CG_TOP` [:strong:`-maxr` :emphasis:`<float>`] [:strong:`-kENM` :emphasis:`<float>`] [:strong:`-pspica`]

Description
-----------
``ENM`` applies the elastic network model (ENM) to CG configuration and generates
topology files to prepare input files for performing CG-MD with SPICA.
Cutoff length and force constant values are required for the ENM and can be 
specified with the ``-maxr`` and ``-kENM`` options, respectively. In default,
They are set to 9.0 A and 1.195 kcal/A2, and these values are employed in the SPICA
protein model (See the `protein SPICA`_ paper).

.. _protein SPICA: https://pubs.acs.org/doi/10.1021/acs.jctc.1c01207

Example
-------
``cg_spica ENM 2mag.cg.pdb 2mag.cg.top`` 

.. image:: figures/ENM_example.jpg
    :scale: 30
    
protein CG configuration (CPK) applied ENM (red)

:download:`2mag.cg.pdb <data/2mag.cg.pdb>` 
:download:`2mag.cg.top <data/2mag.cg.top>`  


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
``-kENM`` <float> (1.195) [kcal/A2]
    force constant of harmonic potential for ENM
``-pspica`` (off)
    Assign partial charge (0.5590) for pSPICA FF (default: 0.1118, for SPICA FF)

