#include "MGDCM.h"
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>
#include <Domain.h>
#include <ID.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <Node.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <TaggedObjectStorage.h>
#include <EquiSolnAlgo.h>

void *OPS_MGDCM()
{
	double lanmbda;
	if (OPS_GetNumRemainingInputArgs() < 1)
	{
		opserr << "WARNING integrator MGDCM needs the initial load factor \n";
		return 0;
	}

	int numdata = 1;
	if (OPS_GetDoubleInput(&numdata, &lanmbda) < 0)
	{
		opserr << "WARNING integrator MGDCM failed to read initial load factor\n";
		return 0;
	}
	if (OPS_GetNumRemainingInputArgs() > 0)
	{
		int maxIter;
		if (OPS_GetIntInput(&numdata, &maxIter) < 0)
		{
			opserr << "WARNING integrator MGDCM failed to read maxIter\n";
			return 0;
		}

		int minIter;
		if (OPS_GetIntInput(&numdata, &minIter) < 0)
		{
			opserr << "WARNING integrator MGDCM failed to read minIter\n";
			return 0;
		}

		int momentum;
		if (OPS_GetIntInput(&numdata, &momentum) < 0)
		{
			opserr << "WARNING integrator MGDCM failed to read momentum\n";
			return 0;
		}
		if (momentum == 1)
		{
			opserr << "Using MGDCM with momentum, maxIter = " << maxIter << ", minIter = " << minIter << endln;
			return new MGDCM(lanmbda, maxIter, minIter, true);
		}
		else
			return new MGDCM(lanmbda, maxIter, minIter, false);
	}
	else
	{
		return new MGDCM(lanmbda, 15, 3, false);
	}
}

MGDCM::MGDCM(double dLambdaBar, int maxIt, int minIt, bool momentum)
	: StaticIntegrator(INTEGRATOR_TAGS_MGDCM),
	  deltaLambdaBar(dLambdaBar), numgsp(0.0),
	  i(0), k(0), dLambda(0.0),
	  dupp1(0), dupc1(0),
	  duHat(0), duBar(0), Fext(0)
{
	this->i = 0;
	this->k = 0;
	this->numgsp = 0.0;
	this->lambda = 0.0;
	this->deltaLambdaBar = dLambdaBar;
	this->maxIterMomentum = maxIt;
	this->minIterMomentum = minIt;
	this->useMomentum = momentum;
	this->sign = 1;
	this->numIterLastStep = 0;
}

void MGDCM::getFext()
{
	AnalysisModel *theModel = this->getAnalysisModel();
	LinearSOE *theLinSOE = this->getLinearSOE();
	if (theModel == 0 || theLinSOE == 0)
	{
		opserr << "WARNING ArcLength::newStep() ";
		opserr << "No AnalysisModel or LinearSOE has been set\n";
		return;
	}
	double currentLambda = theModel->getCurrentDomainTime();
	theModel->applyLoadDomain(1.0);
	if (Fext == 0 || Fext->Size() != theLinSOE->getNumEqn())
	{
		if (Fext != 0)
			delete Fext;
		Fext = new Vector(theLinSOE->getNumEqn());
		if (Fext == 0 || Fext->Size() != theLinSOE->getNumEqn())
		{
			opserr << "MGDCM::getFext - out of memory\n";
			exit(-1);
		}
	}
	theLinSOE->zeroB();
	int res = this->formNodalUnbalance();
	(*Fext) = theLinSOE->getB();
	theModel->applyLoadDomain(currentLambda);
	this->formUnbalance();
}

MGDCM::~MGDCM()
{
	if (duHat != 0)
		delete duHat;
	if (duBar != 0)
		delete duBar;
	if (Fext != 0)
		delete Fext;
}

double MGDCM::getDeltaLambda(Vector dup, Vector dur)
{
	if (k == 1)
	{
		double dl;
		if (i == 1)
		{
			numgsp = (dup ^ dup);
			sign = 1;
			if (numgsp < 0.0)
			{
				sign = -1;
			}
			dl = deltaLambdaBar;
		}
		else
		{
			int nsign = 1;
			float fact = *dupp1 ^ dup;
			if (fact < 0.0)
			{
				nsign = -1;
			}

			sign = nsign * sign;
			double gsp = numgsp / (dup ^ dup);
			dl = sign * deltaLambdaBar * sqrt(gsp);
		}
		dupp1 = new Vector(dup);
		dupc1 = new Vector(dup);
		return dl;
	}
	return -(*dupc1 ^ dur) / (*dupc1 ^ dup);
}

// Called each time analyze command is called
int MGDCM::newStep(void)
{
	i++;
	LinearSOE *theSOE = this->getLinearSOE();
	if (i == 1)
	{
		this->getFext();
	}
	k = 0;
	AnalysisModel *theModel = this->getAnalysisModel();
	if (theModel == 0)
	{
		opserr << "WARNING MGDCM::newStep() ";
		opserr << "No AnalysisModel has been set\n";
		return -1;
	}
	theModel->applyLoadDomain(lambda);

	if (useMomentum && i > 1)
	{ // only from 2nd step onward
		if (numIterLastStep <= minIterMomentum)
		{
			deltaLambdaBar *= 2.0; // double step size
			opserr << "Last step iters: " << numIterLastStep << endln;
			opserr << "Increasing step size to " << deltaLambdaBar << endln;
		}
		else if (numIterLastStep >= maxIterMomentum)
		{
			deltaLambdaBar *= 0.5; // halve step size
			opserr << "Last step iters: " << numIterLastStep << endln;
			opserr << "Decreasing step size to " << deltaLambdaBar << endln;
		}
	}

	numIterLastStep = 0;
	return 0;
}

int MGDCM::update(const Vector &deltaU)
{

	k++;

	AnalysisModel *myModel = this->getAnalysisModel();
	LinearSOE *theSOE = this->getLinearSOE();
	if (myModel == 0 || theSOE == 0)
	{
		opserr << "WARNING MGDCM::update() ";
		opserr << "No AnalysisModel or LinearSOE has been set\n";
		return -1;
	}
	double ld = myModel->getCurrentDomainTime();

	Vector duguino = theSOE->getX();
	if (k == 1)
	{
		duguino = duguino * 0.0;
	}
	this->formTangent();
	theSOE->setB(*Fext);
	int res = theSOE->solve();
	if (res < 0)
	{
		opserr << "MGDCM::update - failed in LinearSOE solve\n";
		return -1;
	}
	if (duHat == 0 || duHat->Size() != theSOE->getNumEqn())
	{
		if (duHat != 0)
			delete duHat;
		duHat = new Vector(theSOE->getNumEqn());
		if (duHat == 0 || duHat->Size() != theSOE->getNumEqn())
		{
			opserr << "MGDCM::update - out of memory\n";
			exit(-1);
		}
	}
	(*duHat) = theSOE->getX();
	double dLambda = this->getDeltaLambda(*duHat, duguino);
	Vector du = *duHat * dLambda + duguino;
	myModel->incrDisp(du);
	if (myModel->updateDomain() < 0)
	{
		opserr << "LoadControl::update - model failed to update for new dU\n";
		return -1;
	}
	lambda += dLambda;
	myModel->applyLoadDomain(lambda);
	theSOE->setX(du);
	numIterLastStep = k;

	return 0;
}

void MGDCM::Print(OPS_Stream &s, int flag)
{
	AnalysisModel *theModel = this->getAnalysisModel();
	if (theModel != 0)
	{
		double currentLambda = theModel->getCurrentDomainTime();
		s << "\t MGDCM - currentLambda: " << currentLambda;
		s << "  deltaLambda: " << deltaLambdaBar << endln;
	}
	else
		s << "\t MGDCM - no associated AnalysisModel\n";
}

int MGDCM::sendSelf(int cTag, Channel &theChannel) { return 0; }
int MGDCM::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) { return 0; }
