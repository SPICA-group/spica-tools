maptraj
===============
Usage   
--------
.. parsed-literal::

    cg_spica maptraj :strong:`AA_PDB` :strong:`AA_TRAJ` :strong:`CG_TRAJ` [:strong:`-json` :emphasis:`<jsonfile>`] [:strong:`-outpdb` <.pdb>] 
    [:strong:`-begin` :emphasis:`<frame>`] [:strong:`-last` :emphasis:`<frame>`] [:strong:`-nodelwat`]

Description
-----------
``maptraj`` executes CG mapping to generate CG trajectory from AA trajectory.

To read and write trajectory files, this program uses `MDAnalysis`_, 
which is a python library to analyze MD simulation trajectories. 
Therefore, available trajectory formats in the program are the same as those 
able to be read with the library. For its installation, see the official web 
page about `MDAnalysis Installation <https://www.mdanalysis.org/pages/installation_quick_start>`_.

The program requires a json formatted file having CG mapping information
that is included in ``spica-tools-devel/src`` directory. 
``-json`` can specify other json files you want to use for mapping.
Option ``-begin`` and ``-last`` can be used to select mapped trajectory frames.

.. _MDAnalysis: https://www.mdanalysis.org


Example
-------
``cg_spica maptraj dopc.aa.pdb dopc.aa.xtc dopc.cg.xtc`` 

Positional args
---------------

``AA_PDB`` [<.pdb>] 
    Input all-atom configuration
``AA_TRAJ`` [<trajectory>] 
    Input all-atom trajectory
``CG_TRAJ`` [<trajectory>] 
    Output coarse-grained trajectory

Optional args
-------------

``-outpdb`` <.pdb> (fin.cg.pdb)
    Output PDB of mapped final CG configuration
``-json`` <jsonfile> (src/spica_top.json)
    Input json file for CG mapping
``-begin`` <int> (0)
    First frame in mapping
``-last`` <int> (-1)
    Last frame in mapping
``-nodelwat`` (off)
    Not delete excess CG water due to CG ion mapping

