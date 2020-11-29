#ifndef _Mesh_
#define _Mesh_





class Mesh
{
public:
	static double Lx, Ly; //size of the domain;
	static int Nx, Ny; //Number of cells;
	double* xc, * yc; //Location of cell centers;
	double* xf, * yf; //Location of cell faces;

	Mesh();
	~Mesh();
};

#endif