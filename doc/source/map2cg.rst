map2cg
===============
Usage   
--------
.. parsed-literal::

    cg_spica map2cg :strong:`AA_PDB` :strong:`CG_PDB` [:strong:`-json` :emphasis:`<jsonfile>`] [:strong:`-nodelwat`] [:strong:`-verbose` :emphasis:`<int>`]

Description
-----------
``map2cg`` performs CG mapping, which converts AA coordinates to CG coordinates.
This program requires a json formatted file with CG mapping information
contained in the ``spica-tools-devel/src`` directory. 
With the ``-json`` option, you can specify other json files to be used for mapping.

Example
-------
``cg_spica map2cg dopc.aa.pdb dopc.cg.pdb`` 

.. image:: figures/map2cg_example.jpg
    :scale: 40

:download:`dopc.aa.pdb <data/dopc.aa.pdb>` 
:download:`dopc.cg.pdb <data/dopc.cg.pdb>`  

Positional args
---------------

``AA_PDB`` [<.pdb>] 
    Input all-atom configuration
``CG_PDB`` [<.pdb>] 
    Output coarse-grained configuration

Optional args
-------------

``-json`` <jsonfile> (src/spica_top.json)
    Input json file for CG mapping
``-nodelwat`` (off)
    Not delete excess CG water due to CG ion mapping
``-verbose`` <int> (1)
    Activate verbose logging (0: off, 1: on)

