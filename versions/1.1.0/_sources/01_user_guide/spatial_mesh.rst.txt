.. _spatial-mesh-section: 

Spatial Mesh 
------------

The spatial mesh can be provided in two different formats. These can be selected with the ``-meshFileType`` keyword in the `control.dat` file, which can be set to ``cartesian`` or ``curvilinear`` for ``.xyz`` and ``.grid`` formats, respectively. While both formats describe a structured mesh (TOSCA can only operate with this type of mesh), the first one is additionally cartesian. Both types of mesh can be defined such that grid points are stretched in all spatial directions. The ``.xyz`` mesh format is TOSCA's custom Cartesian mesh format, and it is very convenient since the user only has to provide the discretization along each of the three cartesian axes. The ``.xyz`` mesh file format can be generalized as

.. code-block:: bash

	Nx Ny Nz
	x(k=1,j=1,i=1) y(k=1,j=1,i=1) z(k=1,j=1,i=1)
	x(k=2,j=1,i=1) y(k=2,j=1,i=1) z(k=2,j=1,i=1)
	:
	x(k=Nk ,j=1,i=1) y(k=Nk ,j=1,i=1) z(k=Nk ,j=1,i=1)
	x(k=1,j=1,i=1) y(k=1,j=1,i=1) z(k=1,j=1,i=1)
	x(k=1,j=1,i=2) y(k=1,j=1,i=2) z(k=1,j=1,i=2)
	:
	x(k=1,j=1,i=Ni) y(k=1,j=1,i=Ni) z(k=1,j=1,i=Ni)
	x(k=1,j=1,i=1) y(k=1,j=1,i=1) z(k=1,j=1,i=1)
	x(k=1,j=2,i=1) y(k=1,j=2,i=1) z(k=1,j=2,i=1)
	:
	x(k=1,j=Nj ,i=1) y(k=1,j=Nj ,i=1) z(k=1,j=Nj ,i=1)

which, for a cube extending 50 m in each direction, with 10 m cells along each direction (hence 6 points), is given by 

.. code-block:: bash

	6 6 6 
	0 0 0
	10 0 0 
	20 0 0
	30 0 0
	40 0 0
	50 0 0
	0 0 0 
	0 10 0 
	0 20 0 
	0 30 0 
	0 40 0
	0 50 0
	0 0 0 
	0 0 10
	0 0 20
	0 0 30
	0 0 40
	0 0 50

As can be noticed, the mesh file is extremely light, even for large meshes, as only the discretization along the three cartesian axes is provided. The coordinates of the remaining mesh cells are inferred iternally within TOSCA, starting from the axes discretization.   

The ``.grid`` format is a curvilinear mesh format (available for example in the meshing software `Pointwise`), where the cartesian coordinates of each mesh point are provided. This becomes useful when the mesh is deformed, and a unique relationship between the cartesian and curvilinear sets of coordinates does not exist anymore. This is a much heavier mesh format, which can be generalized as follows

.. code-block:: bash

	Nj Nk Ni
	x(k=1,j=1,i=1) ..... x(k=1,j=Nj,i=1)
	:
	:
	x(k=Nk,j=1 ,i=1) ....x(k=Nk,j=Nj ,i=1)
	x(k=1,j=1,i=2) ....x(k=1,j=Nj,i=2)
	:
	:
	x(k=Nk ,j=1 ,i=Ni) ...x(k=Nk ,j=Nj ,i=Ni)
	y(k=1,j=1,i=1) ..... y(k=1,j=Nj,i=1)
	:
	:
	y(k=Nk,j=1 ,i=1) ....y(k=Nk,j=Nj ,i=1)
	y(k=1,j=1,i=2) ....y(k=1,j=Nj,i=2)
	:
	:
	y(k=Nk ,j=1 ,i=Ni) ...y(k=Nk ,j=Nj ,i=Ni)
	z(k=1,j=1,i=1) ..... z(k=1,j=Nj,i=1)
	:
	:
	z(k=Nk,j=1 ,i=1) ....z(k=Nk,j=Nj ,i=1)
	z(k=1,j=1,i=2) ....z(k=1,j=Nj,i=2)
	:
	:
	z(k=Nk ,j=1 ,i=Ni) ...z(k=Nk ,j=Nj ,i=Ni)

Taking as an example the same box extending 50 m in each direction, with 10 m cells along each direction, this can be expressed using the ``.grid`` format as

.. code-block:: bash

   6 6 6
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
 
.. note::

   Both in the ``.grid`` and ``.xyz`` format, the user has to provide the *mesh points* rather then the *cell centers* coordinates.
   
Notably, every cartesian mesh can be expressed using the curvilinear format, while the opposite is not true. TOSCA always operates using curvilinear coordinates internally, and since boundary conditions are given on curvilinear patches one has to bear in mind a few concepts, explained below, in order to assign BCs in the way he intends to. Firstly, one can think of curvilinear coordinates as the mesh cell indices in each of the three directions, expressed in TOSCA as `k`, `j` and `i`. With this in mind, it appears clear that each boundary can be identified by keeping one of the indices fixed while the others two vary. In fact, boundaries are identifyed in TOSCA using `kLeft`, `kRight`, `jLeft`, `jRight`, `iLeft` and `iRight`. When the mesh is cartesian (and the axes are aligned to the edges of the domain) boundaries can be also identifyed with the minimum and maximum of each cartesian coordinates. However, this is not general and it would not hold for a curvilinear mesh, or even a cartesian mesh with the domain edges not aligned with the cartesian coordinates. For this reason, boundaries are always identified using the curvilinear definition in TOSCA. To know how to assign boundary conditions, e.g. what boundary corresponds to the `kLeft` patch and so on, one has to keep in mind that the `.grid` file is read by TOSCA using a `i,k,j` nested loop order, namely

.. code-block:: C

   // loop through the mesh cells 
   for (i=0; i<Ni; i++) 
   {
      for (k=0; k<Nk; k++)
      {
         for (j=0; j<Nj; j++)
         {
            // read point k,j,i in the .grid file 
         }
      }
   }
   
Looking at the ``.grid`` file format, it is clear that there is no unique relation between `k,j,i` and `x,y,z`. However, when the mesh is cartesian, this relation depends on how the mesh file is created. It is easy to notice that, in the provided ``.grid`` example, `k,j,i` correspond to `x,z,y`. For the ``.xyz`` file format, where a unique relation exists, TOSCA constructs the mesh such that `k,j,i` correspond to `x,z,y` (this can be also appreciated in the provided example). In an ABL or wind farm simulation, where x is usually aligned with the wind direction, and y and z are the spanwise and vertical coordinates, respectively, the inlet will corrispond to `kLeft`, the outlet to `kRight`, the bottom and top to `jLeft` and `jRight`, respectively, while the spanwise boundaries will be `iLeft` and `iRight`. This relation will always apply as long as a cartesian mesh is used. The relations between curvilinear and cartesian coordinates depending on the mesh and domain type are summarized in the following table. 

+------------------+----------------------------+------------------------------+
|  **format**      | **cartesian domain**       | **curvilinear domain**       |
+------------------+----------------------------+------------------------------+
|   ``.grid``      | case-dependent             | unique relation impossible   |
+------------------+----------------------------+------------------------------+
|   ``.xyz``       | `k,j,i` -> `x,z,y`         | non-representable            |
+------------------+----------------------------+------------------------------+
