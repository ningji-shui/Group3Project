#include "Flow.h"
#include <iostream>

using namespace std;

/* Initialization of basic parameters */
//Scale of the domain [m]
double Mesh::Lx = 2.;
double Mesh::Ly = 1.;

//Number of cells at each dimension
int Mesh::Nx = 128;
int Mesh::Ny = 64;

//Simulation time
double Flow::T = 1.;

//Time step
double Flow::dt = 0.001;

//Reynold number
double Flow::Re = 10000.;

int main()
{
	Flow* flow = new Flow;

		flow->updateFlow();

	return 0;
}