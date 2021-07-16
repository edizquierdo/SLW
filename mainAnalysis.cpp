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
const int N = 2;
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
		Insect.Step1CPG(StepSize);
		pastFootState = newFootState;
		newFootState = Insect.Leg.FootState;
	}

	// Run the agent until you find a crossing point
	double timer = 0.0;
	while (! ((newFootState < 0.5) && (pastFootState > 0.5)))
	{
		Insect.Step1CPG(StepSize);
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
	Insect.Step1CPG(StepSize);
	pastFootState = newFootState;
	newFootState = Insect.Leg.FootState;
	while (! ((newFootState < 0.5) && (pastFootState > 0.5)))
	{
		Insect.Step1CPG(StepSize);
		pastFootState = newFootState;
		newFootState = Insect.Leg.FootState;
		timer += StepSize;
		if (timer > MaxWaitTime)
		{
			return 0.0;
		}
	}
	endingx = Insect.cx;

	// If found, return velocity
	return (endingx-startingx)/timer;
}
// ------------------------------------
// Fitness Maps
// ------------------------------------
int main (int argc, const char* argv[])
{
	// ./2NA 2Ai.best.gen.dat 2Ai.broad.ps.b.csv 3 -16.0 16.0 0.01 4 -16.0 16.0 0.01
	// ./2NA 2Ai.best.gen.dat 2Ai.broad.ps.sc.csv 5 -16.0 16.0 0.01 8 -16.0 16.0 0.01
	// ./2NA 2Ai.best.gen.dat 2Ai.broad.ps.cc.csv 6 -16.0 16.0 0.01 7 -16.0 16.0 0.01

	double fitness;
	ofstream psfile;
	ifstream genefile;
	std::string const & gFileName = argv[1];
	std::string const & psFileName = argv[2];

	double gene1 = atoi(argv[3]);
	double from1 = atof(argv[4]);
	double to1 = atof(argv[5]);
	double div1 = atof(argv[6]);
	double add1 = atof(argv[7]);
	double in1 = atof(argv[8]);

	double gene2 = atoi(argv[9]);
	double from2 = atof(argv[10]);
	double to2 = atof(argv[11]);
	double div2 = atof(argv[12]);
	double add2 = atof(argv[13]);
	double in2 = atof(argv[14]);

	from1 = (from1 - add1)/div1;
	to1 = (to1 - add1)/div1;
	from2 = (from2-add2)/div2;
	to2 = (to2-add2)/div2;

	cout << from1 << " " << to1 << " " << from2 << " " << to2 << endl;

	genefile.open(gFileName);
	psfile.open(psFileName);
	TVector<double> genotype(1, VectSize);
	genefile >> genotype;
	for (double p1 = from1; p1 <= to1; p1 += in1)
	{
		genotype[gene1] = p1;
		for (double p2 = from2; p2 <= to2; p2 += in2)
		{
			genotype[gene2] = p2;
			fitness = AsymptoticFitness(genotype);
			psfile << (p1*div1) + add1 << "," << (p2*div2) + add2  << "," << fitness << endl;
		}
	}
	psfile.close();
	genefile.close();
	return 0;
}
