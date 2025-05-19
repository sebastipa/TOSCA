.. _execution-section: 

Execution
---------

Executing TOSCA is pretty straightforward. Once the code has been compiled (see :ref:`installation_section`), you should see two 
executables in the *build* directory: *tosca* and *tosca2PV*. 

The first executable is the TOSCA solver, while the second allows to convert TOSCA output binary files to *.xmf* files that can be 
read by *ParaView*. 

In order to run TOSCA, copy the two executables to the directory where your input files are located. If running test cases,
we suggest to copy the *tests* directory to a directory of your choiche, in order to maintain your TOSCA repository clean.

`Running in serial` 
~~~~~~~~~~~~~~~~~~~~~

To run in serial, simply run the following command in a terminal:

.. code-block:: bash

    ./tosca 

And to post process the output, run:

.. code-block:: bash

    ./tosca2PV

Finally, open the *xmf* files contained inside the *XMF* directory using *ParaView*.

`Running in parallel` 
~~~~~~~~~~~~~~~~~~~~~~

To run in parallel with **all available** processors (e.g. when submitting jobs on a cluster), simply run the following command in a terminal:

.. code-block:: bash

    mpirun ./tosca 

or with a specific number of processors:

.. code-block:: bash

    mpirun -np 4 ./tosca 

where *4* is the number of processors you want to use. To post process the output in serial (this will also convert 3D fields to *.xmf* format),
run the following command in a terminal:

.. code-block:: bash

    ./tosca2PV

If you want to post process the output in parallel (this will only process the flow sections, usually for large cases where 3D fields are 
only saved for potential restart), run the following command in a terminal:

.. code-block:: bash

    mpirun -np 4 ./tosca2PV 

where *4* is the number of processors you want to use.

Finally, open the *xmf* files contained inside the *XMF* directory with *ParaView*.
