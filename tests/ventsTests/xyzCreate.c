#include <unistd.h>
#include <string.h>
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
    char filename[256]="mesh.xyz";

	// length of the domain in the cartesian x, y and z direction
	double  lengthX = 1.52 , lengthY = 1.52, lengthZ = 1.52;

	// number of nodes in the cartesian x, y and z direction
	int     nx = 76, ny = 76, nz = 76;

	// set periodicity type in the computational i, j and k direction
	int     iPeriodicity = 0, jPeriodicity = 0, kPeriodicity = 0;
	double  dx, dy, dz;

	// reference origin - bottom left co-ordinate of the domain
	vector<double> origin{0, 0, 0};

	dx = lengthX / (nx-1);
	dy = lengthY / (ny-1);
	dz = lengthZ / (nz-1);

	ofstream file;

	file.open(filename);

	if ( file.is_open() )
	{
	    // write the periodicity type if periodic
	    if(iPeriodicity != 0)
	    {
	        file << "-iPeriodicType  " << iPeriodicity << endl;
	    }

	    if(jPeriodicity != 0)
	    {
	        file << "-jPeriodicType  " << jPeriodicity << endl;
	    }

	    if(kPeriodicity != 0)
	    {
	        file << "-kPeriodicType  " << kPeriodicity << endl;
	    }

	    // write the number of POINTS in x y z
	    file << nx << "    "
	         << ny << "    "
	         << nz << "    "
	         << endl;

	    // writing the x axis
	    for(int i=0; i < nx; i++)
	    {
	        file << origin[0]+ i*dx << "    "
	             << 0.0 << "    "
	             << 0.0 << "    "
	             << endl;
	    }

	    // writing the y axis
	    for(int i=0; i < ny; i++)
	    {
	        file << 0.0 << "    "
	             << origin[1]+ i*dy << "    "
	             << 0.0 << "    "
	             << endl;
	    }

	    // writing the z axis
	    for(int i=0; i < nz; i++)
	    {
	        file << 0.0 << "    "
	             << 0.0 << "    "
	             << origin[2]+ i*dz << "    "
	             << endl;
	    }

	    file.close();
	}

	else
	{
	    cout << "Unable to open file";
	}

    return(0);
}
