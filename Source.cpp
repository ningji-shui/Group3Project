#include "Flow.h"
#include <iostream>
#include <fstream>

using namespace std;

/* Initialization of basic parameters */
//Scale of the domain [m]
double Mesh::Lx = 2.;
double Mesh::Ly = 1.;

//Number of cells at each dimension
int Mesh::Nx = 128;
int Mesh::Ny = 64;

//Simulation time
double Flow::T = 1.170;

//Time step
double Flow::dt = 0.001;

//Reynold number
double Flow::Re = 10000.;

int main()
{
	Flow* flow = new Flow;
	flow->updateFlow();
	ofstream fout;
	fout.open("niubility.dat", ios::out | ios::trunc);
	fout << "VARIABLES = \"X\", \"Y\", \"VELOCITY_X\", \"VELOCITY_Y\", \"PRESSURE\"" << endl;
	fout << "ZONE I=" << Mesh::Nx << ", J=" << Mesh::Ny << ", F=POINT" << endl;
	Mesh* mesh1 = new Mesh;
	for (int j = 0; j < Mesh::Ny; j++)
		for (int i = 0; i < Mesh::Nx; i++)
		{
			fout << mesh1->xc[i] << " " << mesh1->yc[j] << " " << 0.5 * (flow->uu[i][j] +
				flow->uu[i + 1][j]) << " " << 0.5 * (flow->vv[i][j] + flow->vv[i][j + 1]) << " " <<
				flow->pp[i][j] << endl;
		}
	fout.close();



	return 0;
}