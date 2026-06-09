import os, sys, time

def readmesh(filename):
    """This will read the Gmsh mesh 'filename' and return
    points and triangles"""
    
    f = open(filename, 'r')

    #################################################################
    #
    # read the first line to get the number of points
    #
    #################################################################

    thisline = f.readline()
    if thisline != "$MeshFormat\n" :
        bail("Unexpected input " + thisline)

    thisline = f.readline()
    segments  = thisline.split( )
    # `.msh' file format, version 2.0
    if segments[0] != "2" and segments[0] != "2.1" and segments[0] != "2.2":
      bail("Unexpected input " + thisline)
    # file-type: ASCII format
    if int(segments[1])!=0 :
      bail("Unexpected input " + thisline)
    # data-size (ignored)

    thisline = f.readline()
    if thisline != "$EndMeshFormat\n" :
        bail("Unexpected input " + thisline)

    #################################################################
    #
    # read the points
    #
    #################################################################

    thisline = f.readline()
    if thisline != "$Nodes\n" :
        bail("Unexpected input " + thisline)
      
    thisline = f.readline()
    number_of_points = int(thisline)

    points = []

    for i in range(number_of_points):
        thisline  = f.readline()
        segments  = thisline.split( )
        thispoint = list(map(lambda p: float(p), segments[1:4]))
        points.append(thispoint)

    thisline = f.readline()
    if thisline != "$EndNodes\n" :
        bail("Unexpected input " + thisline)
      
    #################################################################
    #
    # now read all the volume elements
    # first, ascertain the number of volume elements
    #
    #################################################################

    thisline = f.readline()
    if thisline != "$Elements\n" :
        bail("Unexpected input " + thisline)
      
    thisline             = f.readline()
    number_of_elements   = int(thisline)

    #################################################################
    #
    # then read and populate, keeping track
    # of the tetrahedra subdomain
    #
    #################################################################
    
    triangles           = []
    triangles_subdomain = []

    triangles_id=0
    #tetrahedra_id=0

    for i in range(number_of_elements):
        thisline = f.readline()
        segments = thisline.split( )
        tags = int(segments[2])

        # element type 2: 3 node triangle
        if int(segments[1]) == 2 :
            if len(segments) != 6+tags:
                print ("Found", len(segments), "elements, but expected ", 6+tags)
                print (segments)
                raise Exception("File format error")

            triangles_subdomain.append(triangles_id)
            triangle = list(segments[5:8])
            triangles.append(triangle)
            triangles_id=triangles_id+1

    thisline = f.readline()
    if thisline != "$EndElements\n" :
        bail("Unexpected input " + thisline)
      
    #################################################################
    #
    # send back the read mesh
    #
    #################################################################
    
    return (points, triangles, triangles_subdomain)

def bail(message="Sorry, I don't understand. Bailing out..."):
    """A 'get-out' clause"""
    print (message)
    sys.exit(1)

def initialise_file(file):
    """Attempt to remove the file first to avoid append issues"""
    try:
        os.remove(file)
    except:
        # no problem
        pass

def writeln(line, file):
    """One-shot file line append; this keeps the code terse"""
    f = open(file, 'a')
    f.write(line + '\n')
    f.close()

def stamp():
    """Format a time and a date for stdout feedback"""
    thistime = time.localtime()
    return "[%04d/%02d/%02d %02d:%02d:%02d] " % (thistime[0],
                                                 thistime[1],
                                                 thistime[2],
                                                 thistime[3],
                                                 thistime[4],
                                                 thistime[5])

def writeucd(inputfile, outputfile):
    """Call readmesh to read the Gmsh mesh, then convert and write the UCD mesh"""
    #al = len(sys.argv)
    #if al < 2 or al > 3 or sys.argv[1] == '--help' or sys.argv[1] == '-h':
    #    bail("Usage: " + sys.argv[0] + " gmsh.msh [ucdmesh.inp]")

    #inputfile = sys.argv[1]

    #if len(sys.argv) == 3:
    #    outputfile = sys.argv[2]
    #else:
    #    outputfile = os.path.splitext(inputfile)[0] + ".inp"

    print ()
    print ("About to read, check and process the gmsh file ", inputfile)
    print ("and then create the AVS/UCD mesh", outputfile)
    print ()

    print (stamp() + "Initialising file" + outputfile)
    initialise_file(outputfile)
    print (stamp() + "File " + outputfile + " initialised")

    print (stamp() + "Attemping read of input " + inputfile)
    (points, triangles, trisub) = readmesh(inputfile)
    print (stamp() + "Input file successfully read")

    print (stamp() + "Points:     %d" % len(points))
    print (stamp() + "Triangles:  %d [%d in subdomain]" % (len(triangles),  len(trisub)))

    #################################################################
    #
    # the first line of a UCD file reads:
    # a b c d e
    #
    # where a is the number of nodes
    #       b is the number of cells
    #       c is the length of vector data associated with the nodes
    #       d is the length of vector data associated with the cells
    #       e is the length of vector data associated with the model
    #
    # example: 12 2 1 0 0
    #
    #################################################################

    # n.b. here the third integer indicates the number of placeholders
    
    print (stamp() + "Creating descriptor [UCD]")
    f = open(outputfile, 'w')
    f.write("# UCD geometry file converted from GMSH file " + inputfile + '\n')
    f.write("# File created " + stamp() + '\n')
    f.write("# \n")
    f.write(str(len(points)) + " " + str(len(triangles)) + " 3 0 0" + '\n')
    print (stamp() + "UCD descriptor created")
    
    #################################################################
    #
    # then we have nodes in threespace, one line per node
    # n x y z
    # where n is the node ID -- integer (not necessarily sequential)
    #       x,y,z are the x, y and z coordinates
    #
    #################################################################

    print (stamp() + "Now converting nodes")
    for i in range(len(points)):
        x, y, z = points[i][0], points[i][1], points[i][2]
        f.write(str(i+1) + " " + str(x) + " " + str(y) + " " + str(z) + '\n')
    print (stamp() + "Nodes converted")

    #################################################################
    #
    # now the cells, one line/cell:
    #
    # c m t n1 n2 ... nn
    #
    # where c is the cell ID
    #       m is the material type (int, leave as 1 if we don't care)
    #       t is the cell type (prism|hex|pyr|tet|quad|tri|line|pt)
    #       n1...nn is a node ID list corresponding to cell vertices
    #
    #################################################################

    print (stamp() + "Converting triangles")
    
    cellctr = 0
    for triang in triangles:
        tri = triang[0:3] # just the triangle; triang[1] is the material
        tristring = ""
        for node in tri:
            tristring += " " + str(node)
        f.write(str(cellctr+1) + " " + str(0) + " tri" + tristring + '\n')
        cellctr += 1
    
    print (stamp() + "Triangles converted")
    
    f.close()
    ################################################################

    # all done

    print (stamp() + "All finished")


infile = str(sys.argv[1]) + '.msh'
outfile = str(sys.argv[1]) + '.inp'
writeucd(infile, outfile)