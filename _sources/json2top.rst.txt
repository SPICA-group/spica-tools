json2top
===============
Usage   
--------
.. parsed-literal::

    cg_spica json2top :strong:`RESNAME` [:strong:`-json` :emphasis:`<jsonfile>`] [:strong:`-dupang`] 

Description
-----------
``json2top`` generates the topology files needed to setup CG-MD with SPICA.
The output is ``(RESNAME).top``, where ``RESNAME`` is a mandatory argument to run
the program. 
This command requires a file in json format with CG mapping information,
contained in ``spica-tools/src`` directory. 
``-json`` can specify other json files to be used to generate topology files.
Angle indices in generated topology files are created based on bond information decribed
in the json file. By default, for cyclic molecules, duplicated angle indices are 
automatically removed. 
This removal can be turned off with the ``-dupang`` option. 

**NOTE**: After generating topology files, please make sure that the bond and angle indices 
in the file are correct and as expected, especially for molecules containing cyclic groups.

Example
-------
.. code-block:: bash

    cg_spica json2top POPC

:download:`POPC.top <data/POPC.top>` ::

    atom     1  POPC    NC    NC   87.1647    0.1118  U  
    atom     2  POPC    PH    PH   94.9716   -0.1118  U  
    atom     3  POPC    GL    GL   41.0725    0.0000  U  
    ...                                                 
                                                        
    bond     1     2 # NC PH                            
    bond     2     3 # PH GL                            
    ...                                                 
                                                        
    angle     1     2     3 # NC PH GL                  
    angle     2     3     4 # PH GL EST1                
    ...                                                 
                                                        
Positional args
---------------

``RESNAME`` [<string>] 
    Residue name included in the json file

Optional args
-------------

``-json`` <jsonfile> (src/spica_top.json)
    Input json file 
``-dupang`` (off)
    Not delete duplicated angle indices

