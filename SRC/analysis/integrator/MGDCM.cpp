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
			return new MGDCM(lanmbda, maxIter, minIter, true);
		else
			return new MGDCM(lanmbda, maxIter, minIter, false);
	}
	else
	{
		return new MGDCM(lanmbda, 15, 2, false);
	}
}

MGDCM::MGDCM(double dLambdaBar, int maxIt, int minIt, bool momentum)
	: StaticIntegrator(INTEGRATOR_TAGS_MGDCM),
	  deltaLambdaBar(dLambdaBar), numgsp(0.0), maxIterMomentum(maxIt), minIterMomentum(minIt),
	  i(0), k(0),
	  useMomentum(momentum), dLambda(0.0),
	  dupp1(0), dupc1(0),
	  duHat(0), duBar(0), du(0), Fext(0)
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
	this->formUnbalance();
	(*Fext) = theLinSOE->getB();
	theModel->applyLoadDomain(currentLambda);
}

MGDCM::~MGDCM()
{
	if (duHat != 0)
		delete duHat;
	if (duBar != 0)
		delete duBar;
	if (du != 0)
		delete du;
	if (Fext != 0)
		delete Fext;
}

double MGDCM::getDeltaLambda(Vector dup, Vector duguino)
{
	Vector dur;
	int sign;
	if (k == 1)
	{
		dur = duguino * 0.0;
		if (i == 1)
		{
			numgsp = (dup ^ dup);
			sign = 1;
			if (numgsp < 0.0)
			{
				sign = -1;
			}
			dupp1 = new Vector(dup);
			dupc1 = new Vector(dup);
			return deltaLambdaBar;
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
			dupp1 = new Vector(dup);
			dupc1 = new Vector(dup);
			double dl = sign * deltaLambdaBar * sqrt(gsp);
			return dl;
		}
	}
	else
	{
		dur = duguino;
	}
	return -(*dupc1 ^ dur) / (*dupc1 ^ dup);
}

// Called each time analyze command is called
int MGDCM::newStep(void)
{
	i++;
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
	this->formUnbalance();
	return 0;
}

int MGDCM::update(const Vector &deltaU)
{

	k++;
	opserr << "Entra a update i=" << i << " k=" << k << endln;
	AnalysisModel *myModel = this->getAnalysisModel();
	LinearSOE *theSOE = this->getLinearSOE();
	if (myModel == 0 || theSOE == 0)
	{
		opserr << "WARNING MGDCM::update() ";
		opserr << "No AnalysisModel or LinearSOE has been set\n";
		return -1;
	}
	// En teoria estas dos lineas no son necesarias
	// Aqui el problema es que por alguna razon se está calculando
	// el x casi en Cero. Es como si no se estuviera actualizando el lambda
	// this->formTangent();
	// this->formUnbalance();
	double ld = myModel->getCurrentDomainTime();
	opserr << "lambda antes de update= " << ld << endln;
	Vector duguino = theSOE->getX();
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
	opserr << "duHat " << *duHat << endln;
	double dLambda = this->getDeltaLambda(*duHat, duguino);
	opserr << "dlambda " << dLambda << endln;
	opserr << "duguino " << duguino << endln;

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
