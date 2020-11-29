#include "Ensemble.h"

Ensemble::Ensemble(Flow* flow_)
{
	flow = flow_;
	dt = Flow::dt;
}

Ensemble::~Ensemble()
{

}

void Ensemble::updateEnsemble()
{
	srand(time(0) + rand());
	n_in = m_in * dt;
	particleinsert();
	Runge_Kutta_explicit();
	particlecheck();
}

void Ensemble::particleinsert()
{
	for (int i = 0; i < n_in; i++)
	{
		Particle* temp;
		temp->x = 0.;
		temp->y = 0.45 + 0.1 * rand() / RAND_MAX;
		temp->u = getuleftBC(temp->y);
		temp->v = 0.;

		P.push_back(*temp);
	}
}

void Ensemble::Runge_Kutta_explicit()
{
	int j = 0;
	for (j = 0; j < P.size; j++)
	{
		P[j].x += Runge_Kutta(P[j].u);
		P[j].y += Runge_Kutta(P[j].v);
		P[j].u += Runge_Kutta(1. / St * (getVx(P[j]) - P[j].u));
		P[j].v += Runge_Kutta(1. / St * (getVy(P[j]) - P[j].v));
	}
}

double Ensemble::Runge_Kutta(double q)
{
	double K1 = q;
	double K2 = q + dt / 2 * K1;
	double K3 = q + dt / 2 * K2;
	double K4 = q + dt * K3;

	return (K1 + 2 * K2 + 2 * K3 + K4) * dt / 6.;
}

double Ensemble::getVx(Particle P_)
{

}

double Ensemble::getVy(Particle P_)
{

}

void Ensemble::particlecheck()
{
	int j = 0;
	for (j = 0; j < P.size; j++)
	{
		if (P[j].x > Mesh::Lx || P[j].y > Mesh::Ly || P[j].y < 0. || P[j].x < 0.)
			P.erase(P.begin + j);
	}
}