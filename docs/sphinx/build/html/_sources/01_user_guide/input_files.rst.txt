.. _input-files-section:

Input Files
-----------

This section describes TOSCA's input files and their available entries. Within the code, specific entry types are defined, which are used to define boundary conditions as well as other inputs throughout the code. These are 

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
 

.. toctree::
   :maxdepth: 2
   
   input_files/control.rst
   input_files/boundary.rst
   input_files/ablProperties.rst
   input_files/turbines.rst
   input_files/ibm.rst

..   
   .. include:: input_files/control.rst
   .. include:: input_files/boundary.rst
   .. include:: input_files/ablProperties.rst
   .. include:: input_files/turbines.rst
