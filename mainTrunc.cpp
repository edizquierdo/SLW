#include <iostream>
#include "TSearch.h"
#include "LeggedAgent.h"
#include "CTRNN.h"
#include "random.h"

// Task params
const double StepSize = 0.01;
const double RunTransientDuration = 200.0;
const double Duration = 10000.0;

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

double TruncatedFitnessFunction(TVector<double> &genotype)
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
		// Run the agent
    for (double time = 0; time < RunTransientDuration; time += StepSize) {
        Insect.Step2CPG(StepSize);
    }
		double startingx,endingx;
		startingx = Insect.cx;
    for (double time = 0; time < Duration; time += StepSize) {
        Insect.Step2CPG(StepSize);
    }
    endingx = Insect.cx;
    return (endingx-startingx)/Duration;
}

// ------------------------------------
// Fitness Maps
// ------------------------------------
int main (int argc, const char* argv[])
{
	ifstream genefile("best.gen.dat");
	ofstream outputfile("truncatedfitness.dat");
	TVector<double> genotype(1, VectSize);
	genefile >> genotype;
	outputfile << TruncatedFitnessFunction(genotype) << endl;
	genefile.close();
	outputfile.close();
	return 0;
}
