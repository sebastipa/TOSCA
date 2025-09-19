User Guide
==========

TOSCA uses ASCII input files, organized in files and dictionaries.  
The code provides some level of input checking, meaning that non-recognized inputs are followed by an error message that
lists available possibilities or points the user to possible problems. 

To get a first idea of how to set up and run a TOSCA simulation, the user is advised to first follow the installation instructions
in :ref:`installation_section` and then run one of the provided tutorial cases, as described in :ref:`examples_section`. We tried 
to provide a list of example cases that is comprehensive enough to allow users to always find a tutorial that is close enough to their own needs. 
Notably, all tutorial cases use very coarse/simple meshes, so that they can run on a laptop. The user needs to adapt these cases to their own needs
in order to use TOSCA for production runs. If problems should arise in doing so, please open a `GitHub issue <https://github.com/sebastipa/TOSCA/issues>`_.

TOSCA adopts a standardized case structure. The minimum-required case structure is depicted on the left of the following code-block, 
while the case structure required to run e.g. atmospheric boundary layer (ABL) simulations with potential 
temperature stratification is shown on the right (i.e. with the addition of the ``boundary/T`` and ``ABLProperties.dat`` files).

The principal control file for a TOSCA simulation is the `control.dat` file, located in the case directory (see 
:ref:`control-subsection` for details).

Depending on the type of simulation that one wishes to perform, flags can be activated in the `control.dat`, which prompt TOSCA 
to read additional input files and data. 

.. code-block:: bash
    
   case_name                case_name
   ├── control.dat          ├── control.dat
   ├── mesh.xyz             ├── mesh.xyz
   ├── tosca                ├── ABLProperties.dat
   ├── tosca2PV             ├── tosca
   └── boundary             ├── tosca2PV
       |── U                └── boundary
       └── nut                  ├── U
                                |── nut
                                └── T
    
A complete list of all entries to all TOSCA input files is contained in Sec. :ref:`input-files-section`. Sec. :ref:`spatial-mesh-section`
describes the two mesh formats available in TOSCA and explains how boundary patches are identified using curvilinear coordinates,
in which TOSCA is formulated. Sec. :ref:`overset-section` explains how to set up TOSCA cases with two-way nested domains (i.e. using the 
overset mesh technique) and Sec. :ref:`acquisition-section` details the various outputs that the user can extract from different kinds 
of simulations using TOSCA's acquisition system. Finally, Sec. :ref:`execution-section` shows a typical TOSCA workflow, where the 
solution is first performed using the `tosca` executable, the simulation is then possibily restarted, and `tosca2PV` is used in the end to 
convert output binary data into `.xmf` format to be inspected in *ParaView*. 

.. toctree::
   :maxdepth: 1
   
   01_user_guide/input_files.rst
   01_user_guide/spatial_mesh.rst
   01_user_guide/overset.rst
   01_user_guide/acquisition.rst
   01_user_guide/execution.rst
















 
    
    
    
    
    
    
    
    
    
    
    
    
    
     
