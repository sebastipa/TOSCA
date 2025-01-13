User Guide
==========

TOSCA uses ASCII input files, organized in files and dictionaries.  
The code provides some level of input checking, meaning that non-recognized inputs are followed by an error message that
lists available possibilities. 

TOSCA has a standardized case structure. The minimum-required case structure is depicted on the right of the following figure, 
while the case structure required to run e.g. atmospheric boundary layer (ABL) and wind farm simulations with potential 
temperature stratification is shown on the right (i.e. with the addition of the ``boundary/T`` and ``ABLProperties.dat`` files).
The principal control file for a TOSCA simulation is the `control.dat` file, located in the case directory (see 
:ref:`control-subsection` for details).
Depending on the type of simulation that one wishes to perform, flags can be activated in the `control.dat`, which prompt TOSCA 
to read additional input files and data. These are described in Sec. :ref:`input-files-section`.

.. image:: ../images/structure_1.png
   :align: center
   
Within TOSCA, specific entry types are defined, which are used to define boundary conditions as well as other inputs throughout the code. These are 

1. *bool*       : an integer which is 0 (false) or 1 (true)
2. *integer*    : an integer number
3. *scalar*     : a floating point number (if an integer is provided, this is cast into a floating point number).
4. *vector*     : a group of floating point numbers defined as *(scalar scalar scalar)*, where the scalar **must** be a floating point number. There shouldn't be any space between parantheses and the first and last vector components. 
5. *string*     : a word without spaces. 
6. *dictionary* : a group of additional entries, where the number is dependant on the name of the dictionary, embodied within `{}` parentheses. Dictionaries can containmultiple integer, bool, scalar, vector and string entries. An example is given below. 
           
.. code-block:: C

   dictionary
   {
      boolEntry    1
      integerEntry 10
      scalarEntry  1.0 // or 1
      vectorEntry  (1.0 1.0 1.0)
      stringEntry  someRandomString
   }
   
Inside the ``boundary`` directory, initial and boundary conditions are set in files having the same name as the field they 
describe (see :ref:`boundary-subsection` for details).  

.. include:: user_guide/input_files.rst
.. include:: user_guide/spatial_mesh.rst
















 
    
    
    
    
    
    
    
    
    
    
    
    
    
     
