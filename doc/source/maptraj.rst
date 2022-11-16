maptraj
===============
Usage   
--------
.. parsed-literal::

    cg_spica maptraj :strong:`AA_PDB` :strong:`AA_TRAJ` :strong:`CG_TRAJ` [:strong:`-json` :emphasis:`<jsonfile>`] [:strong:`-outpdb` <.pdb>] 
    [:strong:`-begin` :emphasis:`<frame>`] [:strong:`-last` :emphasis:`<frame>`] [:strong:`-nodelwat`]

Description
-----------
``maptraj`` executes CG mapping to generate CG trajectories from AA trajectories.

This program uses `MDAnalysis`_, a python library for analyzing MD trajectories
to read and write trajectory files. 
Therefore, the file format of MD trajectory available in this program is the same as 
that can be read by this library. For installation instructions, see the official web 
page about `MDAnalysis Installation <https://www.mdanalysis.org/pages/installation_quick_start>`_.

The program requires a json formatted file with CG mapping information
contained in the ``spica-tools/src`` directory. 
The option ``-json`` enables to select other json files to used for mapping.
The option ``-begin`` and ``-last`` can be used to select the frames of the
mapped trajectory.

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

