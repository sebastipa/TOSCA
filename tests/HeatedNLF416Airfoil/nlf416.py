import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial

def LagInterp(X,Y,N,O):
    # Lagrange interpolation sub-routine. 
    # Inputs:
    # X - sorted x values
    # Y - corresponding y values
    # N - number of uniformly spaced values between the input coordinates to interpolate
    # O - order of the Lagrange polynomial interpolant
    # Output is a tuple containing arrays of the interpolated x and y coordinates
    
    X_interp = [] 
    Y_interp = [] 
    for i in range(0, np.size(X)-O+2):   
        poly = lagrange(X[i:i+O+1], Y[i:i+O+1])
        x_temp = np.linspace(X[i], X[i+1],N+2)
        y_temp = Polynomial(poly.coef[::-1])(x_temp)
        X_interp = np.append(X_interp, x_temp[:-1])
        Y_interp = np.append(Y_interp, y_temp[:-1])

    poly = lagrange(X[-O:], Y[-O:])
    x_temp = np.linspace(X[-2], X[-1],N+2)
    y_temp = Polynomial(poly.coef[::-1])(x_temp)
    X_interp = np.append(X_interp, x_temp)
    Y_interp = np.append(Y_interp, y_temp)
    
    return X_interp, Y_interp

def Rotate(X, Y, angle):
    # rotate a vector by angle degrees
    
    alpha = -float(angle)*np.pi/180.0
    Xrot = X*np.cos(alpha) - Y*np.sin(alpha)
    Yrot = X*np.sin(alpha) + Y*np.cos(alpha)
    
    return Xrot, Yrot

def Mag(X, Y):
    # magnitude of a 2D vector given two components
    
    return np.sqrt(X**2 + Y**2)

def Nearest(array, values):
    # return indices corresponding to the nearest locations of a list of values within the array 
    
    index = np.zeros(np.size(values))
    for i in range(0, np.size(index)):
        index[i] = (np.abs(array - values[i])).argmin()
        
    return index.astype(int)

def writeProbeFile(probes, outfile, homoLoc):
    # Write a TOSCA formatted probe locations file
    # Inputs:
    #  probes = array of probe locations
    #  output = output probe file
    #  homoLoc = homogeneous coordinate location (i.e. spanwise position)
    
    f = open(outfile, 'w', newline='\n')
    line = "probesNumber " + str(probes.shape[0]) + "\n"
    f.write(line)
    line = "timeStart 0 \n"
    f.write(line)
    line = "intervalType  timeStep \n"
    f.write(line)
    line = "timeInterval 1.0 \n"
    f.write(line)
    line = "fields U \n"
    f.write(line)
    line = '\n'
    f.write(line)
    line = "locations \n \n"
    f.write(line)
    
    for i in range(probes.shape[0]):
        line = str(probes[i,0]) + " " + str(homoLoc) + " " + str(probes[i,1]) + "\n"
        f.write(line)
    f.close()

def createProbes(coords, surfFlag, chordwise_indices, normal_distance, rakeType):
    # Write a probe file for sampling by TOSCA
    # Inputs:
    #    coords = airfoil coordinates on either PS or SS
    #    surfFlag = surface flag, either 'SS' or 'PS'
    #    chordwise_indices = chordwise index of probe points along the airfoil surface
    #    normal_indices = wall-normal index of probe points away from the airfoil surface
    
    # Define surface tangent vectors
    tangent = np.zeros((np.size(coords[0]),2))
    tangent[:-1,0] = coords[0][1:]-coords[0][:-1]
    tangent[:-1,1] = coords[1][1:]-coords[1][:-1]
    tangent[-1,:] = tangent[-2,:]
    
    # Rotatae 90 degrees to get surface normal vectors
    if surfFlag == 'SS': 
        normal = Rotate(tangent[:,0], tangent[:,1], -90)
    elif surfFlag == 'PS':
        normal = Rotate(tangent[:,0], tangent[:,1], 90)
    else:
        raise Exception("Surface flag must be either SS or PS")
    
    # Normalize to get unit normal vectors
    unit_normal = np.copy(normal)
    unit_normal[0] = normal[0]/Mag(normal[0], normal[1])
    unit_normal[1] = normal[1]/Mag(normal[0], normal[1])    
    
    # Define number of probes in chordwise direction
    num_chordwise_probes = np.size(chordwise_indices)
    
    # Define number of probes in the surface-normal direction
    num_normal_probes = np.size(normal_distance)
    
    # Total number is the product of surface-normal and chordwise probes
    probes = np.zeros((num_chordwise_probes*num_normal_probes,2))
    
    # Define the probe positions by scaling the unit normal vector by the surface-normal distance
    for i in range(0,num_normal_probes):
        l = num_chordwise_probes*i
        ll = num_chordwise_probes*(i+1)
        probes[l:ll,0] = coords[0][chordwise_indices] + unit_normal[0,chordwise_indices]*normal_distance[i]
        probes[l:ll,1] = coords[1][chordwise_indices] + unit_normal[1,chordwise_indices]*normal_distance[i]
    
    # Print individual probe rakes to files
    if rakeType == "xc":
        # Print one probe rake file for each chordwise position
        for i in range(0,num_chordwise_probes):
            outfile = "probes_" + surfFlag + "_xc_eq_" + str(i)
            l = np.arange(i,num_chordwise_probes*num_normal_probes,num_chordwise_probes)
            writeProbeFile(probes[l,:], outfile, 0.2)
    elif rakeType == "zc":
        # Print one probe rake file for each surface-normal position
        for i in range(0,num_normal_probes):
            l = np.arange(i*num_chordwise_probes,(i+1)*num_chordwise_probes,1)
            outfile = "probes_" + surfFlag + "_zc_eq_" + str(i)
            writeProbeFile(probes[l,:], outfile, 0.2)
    else:
        raise Exception("Rake type must be either xc for fixed chordwise position or zc for fixed surface-normal position.")
    
    return probes
    
def writeGeo(SS, PS, outfile, homoSpan, homoCells, homoCoord):
    # Write a GMSH geo file for the airfoil
    # Inputs: 
    #   SS = SS airfoil coordinates 
    #   PS = PS airfoil coordinates
    #   outfile = output gmsh file name
    #   homoSpan = width of homogeneous direction (spanwise)
    #   homoCells = number of cells in the homogeneous direction
    #   homoCoord = coordinate axis of homogeneous direction
       
    f = open(outfile, 'w', newline='\n')
    f.write('//GMSH input file \n')
    f.write('\n')
    f.write('//Suction side coordinates starting from the trailing edge \n')
    
    # SS cell characteristic length
    dx = SS[0][1:]-SS[0][:-1]
    dy = SS[1][1:]-SS[1][:-1]
    l = np.sqrt(np.power(dx,2)+np.power(dy,2))
    l = np.append(l,l[-1])
    
    count = 100
    for i in range(np.size(SS[0])-1,-1,-1):
        line = "Point("
        line += str(count)
        line += ") = { " 
        if homoCoord == 'x':
            line += "0.0, "
            line += np.array2string(SS[0][i]) + ", "
            line += np.array2string(SS[1][i])
        elif homoCoord == 'y':
            line += np.array2string(SS[0][i]) + ", "
            line += "0.0, "
            line += np.array2string(SS[1][i])
        elif homoCoord == 'z':
            line += np.array2string(SS[0][i]) + ", "
            line += np.array2string(SS[1][i]) + ", "
            line += "0.0"
        else:
            print ('Error - homogeneous direction must be a string x, y, or z')
            sys.exit(1)
        line += ', ' + np.array2string(l[i]) + " }; \n"
        f.write(line)
        count += 1
    
    f.write('\n')
    f.write('//Pressure side coordinates starting from the leading edge \n')
    # PS cell characteristic length
    dx = PS[0][1:]-PS[0][:-1]
    dy = PS[1][1:]-PS[1][:-1]
    l = np.sqrt(np.power(dx,2)+np.power(dy,2))
    l = np.append(l,l[-1])
    
    for i in range(0,np.size(PS[0])):
        line = "Point("
        line += str(count)
        line += ") = { " 
        if homoCoord == 'x':
            line += "0.0, "
            line += np.array2string(PS[0][i]) + ", "
            line += np.array2string(PS[1][i])
        elif homoCoord == 'y':
            line += np.array2string(PS[0][i]) + ", "
            line += "0.0, "
            line += np.array2string(PS[1][i])
        elif homoCoord == 'z':
            line += np.array2string(PS[0][i]) + ", "
            line += np.array2string(PS[1][i]) + ", "
            line += "0.0"
        else:
            print ('Error - homogeneous direction must be a string x, y, or z')
            sys.exit(1)
        line += ', ' + np.array2string(l[i]) + " }; \n"
        f.write(line)
        count += 1
    
    line = "Spline(1000) = {100:" + str(count-1) + ", 100};\n"
    f.write(line)
    line = "Line Loop (1) = {1000};\n"
    f.write(line)
    line = "Plane Surface (1) = {1};\n"
    f.write(line)
    if homoCoord == 'x':
        line = 'Extrude {' + str(homoSpan) + ', 0, 0} {Surface{1}; Layers{' + str(homoCells) + '};}\n'
    elif homoCoord == 'y':
        line = 'Extrude {0, ' + str(homoSpan) + ', 0} {Surface{1}; Layers{' + str(homoCells) + '};}\n'
    elif homoCoord == 'z':
       line = 'Extrude {0, 0, ' + str(homoSpan) + '} {Surface{1}; Layers{' + str(homoCells) + '};}\n'
    else:
       print ('Error - homogeneous direction must be a string x, y, or z')
       sys.exit(1)
    f.write(line)
    f.close()    
    print('Done writing to ' + outfile)
    
# Number of additional points between coordinates
N = 4
interpolantOrder = 4
# angle of attack in degrees
AOA = sys.argv[4] 

# Read NLF416 data from  file
data = np.genfromtxt('nlf416.txt', delimiter=",")

# Select suction surface array and sort by increasing values of x
x_SS = data[0:74,0]
y_SS = data[0:74,1]
sorter = np.argsort(x_SS)
X = x_SS[sorter]
Y = y_SS[sorter]

# Interpolate suction surface array
SS_Interp = LagInterp(X, Y, N, interpolantOrder)

# SS chordwise locations of probe points
SS_probe_locs = np.linspace(0.0,1.0,np.floor(1.0/0.02).astype(int)+1)
SS_probe_locs_index = Nearest(SS_Interp[0], SS_probe_locs)

# Rotate 
SS_Interp = Rotate(SS_Interp[0], SS_Interp[1], AOA)

# Select pressure surface array and sort by increasing values of x
x_PS = data[73:,0]
y_PS = data[73:,1]
sorter = np.argsort(x_PS)
X = x_PS[sorter]
Y = y_PS[sorter]

# Interpolate pressure surface array
PS_Interp = LagInterp(x_PS, y_PS, N, interpolantOrder)

# Rotate 
PS_Interp = Rotate(PS_Interp[0], PS_Interp[1], AOA)

# Write interpolated coordinates to file
np.savetxt('nlf416SS.out', SS_Interp, fmt='%.9e', delimiter=' ', newline='\n', header='NLF416 suction surface coordinates. X coordinate is listed first, then Y.', footer='', comments='# ', encoding=None)
np.savetxt('nlf416PS.out', PS_Interp, fmt='%.9e', delimiter=' ', newline='\n', header='NLF416 pressure surface coordinates. X coordinate is listed first, then Y.', footer='', comments='# ', encoding=None)

# Write GMSH geo file
writeGeo(SS_Interp, PS_Interp, "nlf416.geo", sys.argv[1], sys.argv[2], sys.argv[3])

# PS chordwise locations of probe points
#PS_probe_locs = np.linspace(0.02,1.0,np.floor(0.98/0.02).astype(int)+1)
#PS_probe_locs_index = Nearest(PS_Interp[0], PS_probe_locs)

# Wall-normal locations of probe points, relative to surface
#normal_distance = np.array([0.002,0.004,0.006,0.008,0.010,0.012,0.014,0.016,0.018,0.020,0.025,0.03,0.035,0.04,0.05,0.06,0.08,0.1,0.12,0.15,0.2,0.3])

#rakeType = "xc" ## Rake type is "xc" for rake at a fixed chordwise position, or "zc" for rake at a fixed surface-normal position.
#SS_probes = createProbes(SS_Interp, 'SS', SS_probe_locs_index, normal_distance, rakeType)

#normal_distance = np.array([0.002,0.004,0.006,0.008,0.010,0.012,0.014,0.016,0.018,0.020,0.025,0.03,0.035,0.04,0.05,0.06,0.08,0.1,0.12,0.15])
#PS_probes = createProbes(PS_Interp, 'PS', PS_probe_locs_index, normal_distance, rakeType)

# Scatter plot 
#plt.scatter(x_SS, y_SS, label="Orig", s=100)
#plt.scatter(x_PS, y_PS, s=100)
#plt.scatter(SS_Interp[0],SS_Interp[1], label='Suction surface', s=50)
#plt.scatter(PS_Interp[0],PS_Interp[1], label='Pressure surface')
#plt.scatter(SS_probes[:,0],SS_probes[:,1], label = 'SS probes')
#plt.scatter(PS_probes[:,0],PS_probes[:,1], label = 'PS probes')
#plt.axis('square')
#plt.legend()
#plt.show()
