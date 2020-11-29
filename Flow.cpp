#include <math.h>
#include "Flow.h"

Flow::Flow()
{
	dx = Mesh::Lx / Mesh::Nx;
	dy = Mesh::Ly / Mesh::Ny;
	mesh = new Mesh;
	uu = new double* [Mesh::Nx + 1];
	vv = new double* [Mesh::Ny];
	pp = new double* [Mesh::Nx];
	for (int i = 0; i <= Mesh::Nx; i++)
	{
		uu[i] = new double[Mesh::Ny];
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		vv[i] = new double[Mesh::Ny + 1];
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		pp[i] = new double[Mesh::Ny];
	}
	//Initial conditions
	for (int i = 0; i <= Mesh::Nx; i++)
	{
		for (int j = 0; j < Mesh::Ny; j++)
		{
			uu[i][j] = 0.;
		}
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		for (int j = 0; j <= Mesh::Ny; j++)
		{
			vv[i][j] = 0.;
		}
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		for (int j = 0; j < Mesh::Ny; j++)
		{
			pp[i][j] = 0.;
		}
	}
}

Flow::~Flow()
{
	for (int i = 0; i <= Mesh::Nx; i++)
	{
		delete[] uu[i];
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		delete[] vv[i];
	}
	for (int i = 0; i < Mesh::Nx; i++)
	{
		delete[] pp[i];
	}
	delete[] uu;
	delete[] vv;
	delete[] pp;
}

double Flow::u(int i_, int j_)
{
	if (i_ == 0) return exp(-(mesh->yf[j_] - 0.5) * (mesh->yf[j_] - 0.5) / R / R);
	//else if (i_ == -1) return 0.;
	else if (j_ == -1) return uu[i_][j_ + 1];
	else if (j_ == Mesh::Ny) return uu[i_][j_ - 1];
	else if (i_ == Mesh::Nx + 1) return uu[i_ - 1][j_];
	else return uu[i_][j_];
}

double Flow::v(int i_, int j_)
{
	if (i_ == -1) return 0.;
	else if (j_ == 0 || j_ == Mesh::Ny) return 0.;
	else if (i_ == Mesh::Nx) return vv[i_ - 1][j_];
	else return vv[i_][j_];
}

double Flow::p(int i_, int j_)
{
	if (i_ == -1) return pp[i_ + 1][j_];
	else if (j_ == -1) return pp[i_][j_ + 1];
	else if (j_ == Mesh::Ny) return pp[i_][j_ - 1];
	else if (i_ == Mesh::Nx) return pp[i_ - 1][j_];
	else return pp[i_][j_];
}

void Flow::updateFlow()
{
	int i, j;
	for (i = 0; i < Mesh::Nx; i++)
	{
		for (j = 0; j < Mesh::Ny; j++)
		{
			Hx_new = getHx(i, j);
			Hy_new = getHy(i, j);


			Hx_old = Hx_new;
			Hy_old = Hy_new;

		}
	}
}

double Flow::getHx(int i_, int j_)
{
	return
		-u(i_ + 1, j_) - u(i_ - 1, j_) * (u(i_ + 1, j_) + u(i_ - 1, j_) + 2. * u(i_, j_)) / (4. * dx)
			+ (u(i_, j_) + u(i_, j_ - 1) * (v(i_ - 1, j_) + v(i_, j_))
				- (u(i_, j_ + 1) + u(i_, j_)) * (v(i_ - 1, j_ + 1) + v(i_, j_ + 1))) / (4. * dy);
}

double Flow::getHy(int i_, int j_)
{
		
}