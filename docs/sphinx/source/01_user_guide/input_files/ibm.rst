.. _ibm-section:

`IBM`
~~~~~

This section describes how the immersed boundary method (IBM) can be activated in TOSCA. The code supports both static or dynamic 
IBM, where the latter means that the body can move without limitations inside the computational domain. The IBM has been fully parallelized 
and it features fast algorithms to perform both topology tests (inside-outside) or the search that finds, at every iteration, the 
closest fluid mesh cell for each IBM surface triangle. TOSCA's IBM is referred to as a sharp-interface method, where boundary conditions 
are not applied my means of momentum source terms, but rather by directly recunstructing the boundary condition at the IBM fluid cells. 
These cells are selected by using the following algorithm. First, all cells whose center lies inside the IBM body are marked as IBM solid cells.
These cells are blanked when solving the governing equations. Then, wherever an IBM solid cell is adjacent to a fluid cell, this fluid cell is 
transformed into an IBM fluid cell, where BCs are reconstructed, so that IBM solid cells are always enclosed in a layer of IBM fluid cells. 
This distinction, which can be seen as a cell coloring algorithm, is saved in the *nv* field, which marks each cell based on its IBM 
coloring. For instance, *nv* is 0, 1 and 2 for fluid cells, IBM fluid cells and IBM solid cells, respectively. 

Both boundary conditions and wall models can be applied to the IBM fluid cells. BCs are reconstructed by first creating a fictitious 
ghost point on the surface normal (pointing outward) through the center of the fluid cell. The field of interest is then interpolated
here using trilinear interpolation from the surrounding cells. This is then used to apply the boundary condition at the fluid cell, also
considering the velocity of the IBM body. For wall models, teoretical laws of the wall are used, while for genral BCs a linear interpolation 
is employed. Optionally, one can also use two fictitious ghost points, which lead to a quadratic interpolation instead (less stable).

The body surface mesh is read from an UCD file, which is a standard format for surface meshes (this format can be generated e.g. by PointWise). 
Notably, to read the UCD file in ParaView, the extension must me changed to *.inp*. Finally, in TOSCA only the name of the file should be specified,
thereby dropping the extension. Below is an example of a UCD file for a cuboid body, where each fase is split in two triangles: 

.. code-block:: C

    # UCD Geometry file crated by custom python script 
    #
    #
    8 12 0 0 0
    1 -6.83333e+02 -4.33333e+02 0.00000e+00
    2  6.83333e+02 -4.33333e+02 0.00000e+00
    3  6.83333e+02  1.13333e+03 0.00000e+00
    4 -6.83333e+02  1.13333e+03 0.00000e+00
    5 -6.83333e+02 -4.33333e+02 1.20000e+02
    6  6.83333e+02 -4.33333e+02 1.20000e+02
    7  6.83333e+02  1.13333e+03 1.20000e+02
    8 -6.83333e+02  1.13333e+03 1.20000e+02
    1  0 tri 1 2 3
    2  0 tri 1 3 4
    3  0 tri 5 6 7
    4  0 tri 5 7 8
    5  0 tri 1 2 6
    6  0 tri 1 6 5
    7  0 tri 2 3 7
    8  0 tri 2 7 6
    9  0 tri 3 4 8
    10 0 tri 3 8 7
    11 0 tri 4 1 5
    12 0 tri 4 5 8

Notably, there are 8 vertices and 12 triangles in the example above. 

In order to activate the IBM in TOSCA, the ``-ibm`` flag has to be set to 1 in the `control.dat`` file. TOSCA is prompted to read additional 
inputs inside the ``IBM`` directory. These are the `IBMProperties.dat` file and at least one IBM body in the above format, whose name is 
specified in the `IBMProperties.dat` file. 

IBMProperties.dat 
*****************

The following tables summarize all available entries for each of the `IBMProperties.dat` file. 

.. table:: 
   :widths: 30, 20, 50
   :align: center
                                                                                                       
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**   
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``NumberofBodies``             integer             Number of bodies.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``debug``                      bool                If set to 1, the code will print additional information about the IBM. 
                                                      This is useful for debugging purposes.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``dynamic``                    bool                0 for static IBM, 1 for dynamic IBM.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``checkNormal``                bool                if set to 1 checks that all normals are pointing outward and flips them if 
                                                      necessary. This option is always recommended. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``averageNormal``              bool                averages the normal across more IBM triangles, recommended for irregular  
                                                      bodies or bodies that are under-resolved by the triangular mesh. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``computeForce``               bool                computes forces on the IBM body. Moments are only calculated if the IBM is 
                                                      dynamic. This will be corrected in the future. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``writeSTL``                   bool                writies the IBM body in STL format based on the `writeSettings`. 
                                                      Useful with dynamic IBM, to visualize the body motion.    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``wallShear``                  bool                directly computes the wall shear and applies within the momentum equation.
                                                      Recommended for terrains. Requires `interpolationDistance`. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``interpolationDistance``      scalar              Distance from the body surface where the ghost point is placed. This 
                                                      parameter should be set according to the mean mesh size at the wall and it 
                                                      it is uniform to avoid spikes in the stress produced by variable IBM fluid 
                                                      cell to surface distance. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``abl``                        bool                if set to 1 allows to set `groundLevel`. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``groundLevel``                scalar              Allows to set a user-defined ground level. This used at differet levels in 
                                                      the code when evaluating heigh-dipendend BCs and models. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``writeSettings``              dictionary          Required when `computeForce` is set to 1. Contains entries for the IBM force 
                                                      writing. It is defined as follows: 
                                                      
                                                      .. code-block:: C
                                                    
                                                         writeSettings
                                                         {
                                                            timeStart           scalar // start time for writing
                                                            intervalType        string // timeStep or adjustableTime
                                                            timeInterval        scalar // number of iters or seconds
                                                         }
   ------------------------------ ------------------- ---------------------------------------------------------------------------- 
   ``InterpolationMethod``        string              `CURVIB` (preferred, requires `CURVIBInterpolationType`) or `MSL` 
                                                      (Moving Least Squares). The first method is preferred while the last is 
                                                      only recommended for simple geometries and for no-slip BCs only. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``CURVIBInterpolationType``    string              `CurvibTrilinear` (recommended for terrain) or `CurvibTriangular`. The first 
                                                      interpolates the velocity at the ghost node using tri-linera interpolation 
                                                      from the 8 surrounding fluid cells. The second uses a triangular 
                                                      interpolation from the 3 surrounding fluid cells. The first method is 
                                                      preferred, although it is a bit more expensive. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``interpolationOrder``         string              Only required for `CurvibTrilinear`, can be `linear` (more stable) or 
                                                      `quadratic` (less stable). The latter is recommended for smooth terrain or 
                                                      smooth static bodies. It referst to how many fictitious ghost points are 
                                                      projected from the surface mesh into the fluid domain. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``object<N>``                  dictionary          specifies the object-specific properties. There are as meny dictionaries as 
                                                      the number of bodies. The name of the dictionary is `object` followed by the 
                                                      body index. The body index is an integer starting from 0.
   ============================== =================== ============================================================================
                   
Object Description 
******************

The object-specific properties are defined in the ``object<N>`` dictionary. The following table summarizes the mandatory 
entries that are required for each object when the `bodyMotion` flag is set to *static*. For non-static body motion, 
additional entries are required, which are summarized in the table below.

Static IBM 
----------

.. table:: 
   :widths: 30, 20, 50
   :align: center
                                                                                                       
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**   
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``bodyName``                   string              Name of the body. This is the name of the file without the extension. 
                                                      The file must be in the `IBM` directory. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``bodyType``                   string              This is the type of body, can be `closedBody` or `surfaceBody`. In the 
                                                      first case, the body is a closed watertight body, in the second case, the 
                                                      is one or more open surfaces, and additional entries are required.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``numSurfaces``                integer             Number of surfaces that make up the body. Only required for `surfaceBody`.  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``fileType``                   string or           Can be * ucd*, *grd*, *ascii*, *inp* or *ucd2* for `closedBody`. Can be   
                                  vector of string    * ucd*, *grd*, *inp* or *ucd2* for `surfaceBody`. For `closedBody`, this is 
                                                      a single string, while for `surfaceBody` this is a vector of strings, .e.g 
                                                      (*ucd* *grd* *inp*). The size of the vector must be equal to `numSurfaces`.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``surfaceNames``               vector of string    Names of the surface files that make up a body. This is a vecotor whose 
                                                      elements are the names of the files. The size of the vector must be equal to
                                                      `numSurfaces`. The files must be in the `IBM` directory.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``elementSet``                 string              Used for both `closedBody` and `surfaceBody` when the ``fileType`` is set to 
                                                      *inp* to select an element set from the file. Currently is unique for all 
                                                      surfaces when there are multiple surfaces. For example, given a finite-
                                                      element model mesh, with strings and longerons, it can be used to only 
                                                      select the outer surface of the body for the fluid calculation.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``baseLocation``               vector              Translation vector that allows to translate each point coordinate by this 
                                                      amount in the x, y and z direction. This is useful to move the body directly 
                                                      within TOSCA if this is defined with a different coordinate system.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``bodyMotion``                 string              Type of body motion, can be `static`, `rotation`, `sinusoidal` or
                                                      `pitchingOscillation`. Moving cases required additional entries (see table 
                                                      below). 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``searchCellRatio``            scalar              Factor that multiplies the average cell size, indicating how far away from a 
                                                      given surface triangle to go, when performing the search for the closest 
                                                      fluid cell. A value of 3 is recommended for most cases.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``velocityBCSetType``          string              mode used to set the velocity boundary condition. 
                                                      
                                                      - `matchUiLeft`: match iLeft patch boundary condition 
                                                      - `matchUiRight`: match iRight patch boundary condition 
                                                      - `matchUjLeft`: match jLeft patch boundary condition 
                                                      - `matchUjRight`: match jRight patch boundary condition
                                                      - `setHere`: set in the *IBMProperties.dat* file, requires 
                                                        additional entries. 

                                                      Match-type boundary conditions are useful for terrain simulations, where the 
                                                      IBM body may merge with the bottom patch, in order to ensure that the same 
                                                      boundary condition is applied. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``temperatureBCSetType``       string              mode used to set the temperature boundary condition. Requires `-potentialT`
                                                      to be set to 1 in the `control.dat` file.
                                                      
                                                      - `matchTiLeft`: match iLeft patch boundary condition 
                                                      - `matchTiRight`: match iRight patch boundary condition 
                                                      - `matchTjLeft`: match jLeft patch boundary condition 
                                                      - `matchTjRight`: match jRight patch boundary condition
                                                      - `setHere`: set in the *IBMProperties.dat* file, requires 
                                                        additional entries. 

                                                      Match-type boundary conditions are useful for terrain simulations, where the 
                                                      IBM body may merge with the bottom patch, in order to ensure that the same 
                                                      boundary condition is applied.  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``velocityBC``                 string              Only required with `velocityBCSetType` set to *setHere*.
                                                      Name of the IBM velocity boundary condition. Available entries are 
                                                      `noSlip`, `slip` and `velocityWallFunction`. The latter requires additional 
                                                      entries, depending on `wallFunctionTypeU`, detailed in the next table. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
    ``wallFunctionTypeU``          string             Requires `velocityBC` to be set to `velocityWallFunction`.
                                                      Used to select the type of wall function for the velocity 
                                                      (-1, -3, -4 or -5). 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``temperatureBC``              string              Only required with `temperatureBCSetType` set to *setHere*.
                                                      Name of the IBM velocity boundary condition. Available entries are 
                                                      `zeroGradient`, `fixedValue` (requires `fixedValueT` entry) and 
                                                      `thetaWallFunction`. The latter requires additional entries, depending on 
                                                      `wallFunctionTypeT`, detailed in the next table.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``wallFunctionTypeT``          string              Requires `temperatureBC` to be set to `thetaWallFunction`.
                                                      Used to select the type of wall function for the potential temperature 
                                                      (-2, -3 or -4). 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
    ``fixedValueT``                scalar             Fixed value for the temperature BC. Only required with 
                                                      `temperatureBCSetType` set to *setHere* and `temperatureBC` set 
                                                      to `fixedValue`.
   ============================== =================== ============================================================================

The following table summarizes the velocity wall functions that can be selected from the *IBMProperties.dat* file, together with their 
respective entries. 

.. table:: 
   :widths: 30, 20, 50
   :align: center
                                                                                                       
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**   
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *wallFunctionTypeU -1* (Cabot wall model)
   -------------------------------------------------------------------------------------------------------------------------------    
   ``roughness``                  scalar              equivalent roughness height in meters. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kappa``                      scalar              von Karman constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
    *wallFunctionTypeU -3* (Shumann wall model)
   -------------------------------------------------------------------------------------------------------------------------------    
   ``uStarEval``                  string              `averaged` for laterally-homogeneous flows, `localized` otherwise. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kappa``                      scalar              von Karman constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``thetaRef``                   scalar              reference potential temperature in Kelvin.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``roughness``                  scalar              equivalent roughness height in meters.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gammaM``                     scalar              Shumann model constant.   
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *wallFunctionTypeU -4* (Power law wall model)
   -------------------------------------------------------------------------------------------------------------------------------    
   ``roughness``                  scalar              equivalent roughness height in meters.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kappa``                      scalar              von Karman constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *wallFunctionTypeU -5* (Log law wall model)
   -------------------------------------------------------------------------------------------------------------------------------    
   ``roughness``                  scalar              eqivalent roughness height in meters.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kappa``                      scalar              von Karman constant.
   ============================== =================== ============================================================================

The following table summarizes the temperature wall functions that can be selected from the *IBMProperties.dat* file, together with their 
respective entries. 

.. table:: 
   :widths: 30, 20, 50
   :align: center
                                                                                                       
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**   
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *wallFunctionTypeT -2* (Shumann wall model - wall heat flux prescribed)
   -------------------------------------------------------------------------------------------------------------------------------    
   ``qWall``                      scalar              wall heat flux in J/m2.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
    *wallFunctionTypeT -3* (Shumann wall model - constant heating rate)
   -------------------------------------------------------------------------------------------------------------------------------    
   ``uStarEval``                  string              `averaged` for laterally-homogeneous flows, `localized` otherwise. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kappa``                      scalar              von Karman constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``thetaRef``                   scalar              reference potential temperature in Kelvin.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kRough``                     scalar              equivalent roughness height in meters.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gammaM``                     scalar              Shumann model constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gammaH``                     scalar              Shumann model constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``alphaH``                     scalar              Shumann model constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``heatingRate``                scalar              heating rate in K/s.  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *wallFunctionTypeT -4* (Shumann wall model - variable heating rate)
   -------------------------------------------------------------------------------------------------------------------------------    
   ``uStarEval``                  string              `averaged` for laterally-homogeneous flows, `localized` otherwise. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kappa``                      scalar              von Karman constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``thetaRef``                   scalar              reference potential temperature in Kelvin.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``kRough``                     scalar              equivalent roughness height in meters.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gammaM``                     scalar              Shumann model constant.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gammaH``                     scalar              Shumann model constant.            
   ============================== =================== ============================================================================

Theta wall function type -4 requires look up tables of time (s), surface temperature (K), and Obhukhov length (m), 
stored in ``inflowDatabase/mesoscaleData/time``, ``inflowDatabase/mesoscaleData/surfTemp`` and ``inflowDatabase/mesoscaleData/L``, respectively. 
All vectors should have the same size.

Non-Static IBM 
--------------

The first two entries in the table are required, for each object, for all kinds of non-static body motion. The remaining entries are 
specific to the type of motion listed in the header. 

.. table:: 
   :widths: 30, 20, 50
   :align: center
                                                                                                       
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**   
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``procBoundCenter``            vector              center of the box enclosing the body and its motion throughout the entire
                                                      simulation. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``procBoundSize``              vector              size of the box enclosing the body and its motion throughout the entire
                                                      simulation.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *rotation*
   -------------------------------------------------------------------------------------------------------------------------------    
   ``angularSpeed``               scalar              angular speed in rpm.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``angularAcceleration``        scalar              angular acceleration in rpm per second. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rotationAxis``               vector              rotation axis (normalized by TOSCA).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rotationCenter``             vector              base of the rotation axis. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``maxTipRadius``               scalar              maximum radius of the body, used to compute an equivalent tip speed far 
                                                      from the center of rotation, thus limiting the time step so that the 
                                                      fastest triangular elements of the body only crosse one cell per time step. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
    *sinusoidal*
   ------------------------------------------------------------------------------------------------------------------------------- 
   ``amplitude``                  scalar              motion amplitude in meters. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``frequency``                  scalar              frequency in Hz.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``motionDirection``            vector              direction of the motion (normalized by TOSCA).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *pitchingOscillation*
   -------------------------------------------------------------------------------------------------------------------------------    
   ``angularAmplitude``           scalar              angular amplitude in degrees.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``frequency``                  scalar              frequency in Hz.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``initialAngle``               scalar              initial angle in degrees.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchingAxis``               vector              pitching axis (normalized by TOSCA).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchingCenter``             vector              base of the pitching axis.
   ============================== =================== ============================================================================
