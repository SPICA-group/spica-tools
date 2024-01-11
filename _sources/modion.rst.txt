modion      
===============
Usage   
--------
.. parsed-literal::

    cg_spica modion :strong:`INPUT_PDB` :strong:`OUTPUT_PDB` [:strong:`-conc` :emphasis:`<float>`] [:strong:`-soluq` :emphasis:`<int>`] [:strong:`-pspica`]

Description
-----------
``modion`` adjusts CG NaCl concentration in simulation systems for SPICA or pSPICA.
This program recognizes atom names **W** and **W1** as SPICA and pSPICA CG water models, respectively.
For ions, the atom names **SOD** and **CLA** are recognized as SPICA CG ions,
while the atom names **SOD1** and **CLA1** are recognized as pSPICA CG ions.
The atom names are given in a CG configuration file in the PDB format.

The ``-conc`` option takes a value of the adjusted NaCl concentration after conversion.
The ``-soluq`` option takes a total charge of solute molecules contained in the simulation system.
Note that, in SPICA and pSPICA topology files, the point charges, q = +-1, are set to q = +-0.1118 and q = +-0.5590, 
respectively, to account for the dielectric constant applied in the models.
For an argument of ``-soluq``, specify a total charge of solute molecules with the point charges of q = +-1.
The argument corresponds to the number of counterions needed for neutralization.
The ``-pspica`` option is needed when this program is applied for pSPICA.

Example
-------
.. code-block:: bash

    cg_spica modion prot_memb.cg.pdb prot_memb.cg.200mM.pdb -soluq 4 -conc 0.2

.. code-block:: bash

    cg_spica modion prot_memb.cg.200mM.pdb prot_memb.cg.500mM.pdb -soluq 4 -conc 0.5

.. image:: figures/modion_fig.jpg
    :scale: 15

:download:`prot_memb.cg.pdb <data/setup_gmx_prot/prot_memb.cg.pdb>` 

:download:`prot_memb.cg.200mM.pdb <data/setup_gmx_prot/prot_memb.cg.200mM.pdb>` 
:download:`prot_memb.cg.500mM.pdb <data/setup_gmx_prot/prot_memb.cg.500mM.pdb>` 

Positional args
---------------

``INPUT_PDB`` [<.pdb>] 
    Input configuration file in PDB format
``OUTPUT_PDB`` [<.pdb>] 
    Output configuration file in PDB format

Optional args
---------------

``-conc`` <float> (0.15) [M]
    Specify NaCl concentration to be reached
``-soluq`` <int> (0) [e]
    A total charge of solute molecules (the number of counterions for neutralization)
``-pspica`` (off)
    Apply this program for pSPICA water and ion models (default: for SPICA water and ion models)

