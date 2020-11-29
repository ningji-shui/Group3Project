#include "Mesh.h"

Mesh::Mesh()
{
	//Allocate memory for cell centers/faces locations
	double dx = Lx / Nx;
	double dy = Ly / Ny;
	xf = new double[Nx + 1];
	yf = new double[Ny + 1];
	for (int i = 0; i <= Nx; i++)
	{
		xf[i] = i * dx;
	}
	for (int i = 0; i <= Ny; i++)
	{
		yf[i] = i * dy;
	}

	xc = new double[Nx];
	yc = new double[Ny];
	for (int i = 0; i < Nx; i++)
	{
		xc[i] = xf[i] + dx / 2.;
	}
	for (int i = 0; i < Ny; i++)
	{
		yc[i] = yf[i] + dy / 2.;
	}
}

Mesh::~Mesh()
{
	delete[] xf;
	delete[] yf;
	delete[] xc;
	delete[] yc;
}