json2top
===============
Usage   
--------
.. parsed-literal::

    cg_spica json2top :strong:`RESNAME` [:strong:`-json` :emphasis:`<jsonfile>`] [:strong:`-dupang`] 

Description
-----------
``json2top`` generates topology files needed to setup CG-MD with SPICA.
Output will be ``(RESNAME).top``, where ``RESNAME`` is a mandatory argument to excute
the program. 
This command requires a json formatted file having CG mapping information
that is included in ``spica-tools-devel/src`` directory. 
``-json`` can specify other json files you want to use for generating topology files.
Angle indices in generated topology files are created based on bond information decribed
in the json file. In default, the program automatically deletes duplicated angle indices 
in case of ring molecules. This deletion can be turned off with the ``-dupang`` option. 

**NOTE**: After generating topology files, please make sure the bond and angle indices 
in the file are correct and what you expect, especially for molecules including ring groups.

Example
-------
``cg_spica json2top POPC`` 

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

