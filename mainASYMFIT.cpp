#include <iostream>
#include "TSearch.h"
#include "LeggedAgent.h"
#include "CTRNN.h"
#include "random.h"

// Task params
const double StepSize = 0.01;
const double RunTransientDuration = 200.0;
const double MaxWaitTime = 500.0;

// Nervous system params
const int N = 3;
const double WR = 16.0;
const double BR = 16.0;
const double TMIN = 0.5;
const double TMAX = 10.0;

int	VectSize = N*N + 2*N;

// ------------------------------------
// Genotype-Phenotype Mapping Functions
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
	int k = 1;
	// Time-constants
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), TMIN, TMAX);
		k++;
	}
	// Bias
	for (int i = 1; i <= N; i++) {
		phen(k) = MapSearchParameter(gen(k), -BR, BR);
		k++;
	}
	// Weights
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			phen(k) = MapSearchParameter(gen(k), -WR, WR);
			k++;
		}
	}
}

// ------------------------------------
// Fitness function
// ------------------------------------
double AsymptoticFitness(TVector<double> &genotype)
{
	// Map genootype to phenotype
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	GenPhenMapping(genotype, phenotype);

	// Create the agent
	LeggedAgent Insect;

	// Instantiate the nervous system
	Insect.NervousSystem.SetCircuitSize(N);
	int k = 1;
	// Time-constants
	for (int i = 1; i <= N; i++) {
		Insect.NervousSystem.SetNeuronTimeConstant(i,phenotype(k));
		k++;
	}
	// Bias
	for (int i = 1; i <= N; i++) {
		Insect.NervousSystem.SetNeuronBias(i,phenotype(k));
		k++;
	}
	// Weights
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			Insect.NervousSystem.SetConnectionWeight(i,j,phenotype(k));
			k++;
		}
	}

	Insect.Reset(0, 0, 0);

	// Run the agent for a transient period
	double pastFootState,newFootState;
	for (double time = 0; time < RunTransientDuration; time += StepSize) {
		Insect.StepCPG(StepSize);
		pastFootState = newFootState;
		newFootState = Insect.Leg.FootState;
	}

	// Run the agent until you find a crossing point
	double timer = 0.0;
	while (! ((newFootState < 0.5) && (pastFootState > 0.5)))
	{
		Insect.StepCPG(StepSize);
		pastFootState = newFootState;
		newFootState = Insect.Leg.FootState;
		timer += StepSize;
		if (timer > MaxWaitTime)
		{
			return 0.0;
		}
	}

	// Run the agent until you find the crossing point again
	timer = 0.0;
	double startingx,endingx;
	startingx = Insect.cx;
	Insect.StepCPG(StepSize);
	pastFootState = newFootState;
	newFootState = Insect.Leg.FootState;
	ofstream outputfile("beh.csv");
	outputfile << Insect.NervousSystem.NeuronOutput(1) << "," << Insect.NervousSystem.NeuronOutput(2) << "," << Insect.NervousSystem.NeuronOutput(3) << endl;
	while (! ((newFootState < 0.5) && (pastFootState > 0.5)))
	{
		Insect.StepCPG(StepSize);
		outputfile << Insect.NervousSystem.NeuronOutput(1) << "," << Insect.NervousSystem.NeuronOutput(2) << "," << Insect.NervousSystem.NeuronOutput(3) << endl;
		pastFootState = newFootState;
		newFootState = Insect.Leg.FootState;
		timer += StepSize;
		if (timer > MaxWaitTime)
		{
			return 0.0;
		}
	}
	outputfile.close();
	endingx = Insect.cx;

	// If found, return velocity
	cout << endingx-startingx << endl;
	return (endingx-startingx)/timer;
}
// ------------------------------------
// Fitness Maps
// ------------------------------------
int main (int argc, const char* argv[])
{
	ifstream genefile;
	genefile.open("best.gen.dat");
	TVector<double> genotype(1, VectSize);
	genefile >> genotype;
	AsymptoticFitness(genotype);
	genefile.close();
	return 0;
}
