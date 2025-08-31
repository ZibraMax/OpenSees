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

void *OPS_ArcLength()
{
	double arcLength;
	if (OPS_GetNumRemainingInputArgs() < 2)
	{
		opserr << "WARNING integrator ArcLength arcLength alpha \n";
		return 0;
	}

	int numdata = 1;
	if (OPS_GetDoubleInput(&numdata, &arcLength) < 0)
	{
		opserr << "WARNING integrator ArcLength failed to read arc length\n";
		return 0;
	}
	if (OPS_GetNumRemainingInputArgs() > 0)
	{
		double alpha;
		if (OPS_GetDoubleInput(&numdata, &alpha) < 0)
		{
			opserr << "WARNING integrator ArcLength failed to read alpha\n";
			return 0;
		}
		return new MGDCM(arcLength, alpha);
	}
	else
	{
		return new MGDCM(arcLength);
	}
}

MGDCM::MGDCM(double dLambdaBar, int maxIt, int minIt, bool momentum)
	: StaticIntegrator(INTEGRATOR_TAGS_MGDCM),
	  deltaLambdaBar(dLambdaBar),
	  maxIterMomentum(maxIt), minIterMomentum(minIt),
	  useMomentum(momentum),
	  currentLambda(0.0), dLambda(0.0), signLastStep(1.0),
	  duHat(0), duBar(0), du(0), phat(0)
{
}

int MGDCM::newStep(void)
{
	// assemble load vector for unit lambda
	AnalysisModel *theModel = this->getAnalysisModel();
	LinearSOE *theLinSOE = this->getLinearSOE();
	theModel->applyLoadDomain(1.0);
	phat = &(theModel->getRHS());

	// solve for duHat = K^-1 * phat
	theSOE->setB(*phat);
	if (theSOE->solve() < 0)
		return -1;
	duHat = &(theSOE->getX());

	// initial guess for dLambda
	dLambda = deltaLambdaBar;

	// displacement increment
	if (du == 0)
		du = new Vector(duHat->Size());
	(*du) = (*duHat) * dLambda;

	// apply to model
	theModel->incrDisp(*du);
	return 0;
}

int MGDCM::update(const Vector &deltaU)
{
	// solve for correction direction
	theSOE->setB(deltaU);
	if (theSOE->solve() < 0)
		return -1;
	duBar = &(theSOE->getX());

	// compute new dLambda using MGDCM formula
	double num = (*duHat) ^ deltaU;
	double den = (*duHat) ^ (*duHat);
	dLambda = num / den;

	// total update
	(*du) = (*duHat) * dLambda + (*duBar);

	// apply to model
	theModel->incrDisp(*du);
	currentLambda += dLambda;

	return 0;
}

int MGDCM::domainChanged(void)
{
	// resize vectors to system size
	int sz = theModel->getNumEqn();
	if (du)
		delete du;
	du = new Vector(sz);
	return 0;
}

int MGDCM::sendSelf(int cTag, Channel &theChannel) { return 0; }
int MGDCM::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) { return 0; }

void MGDCM::setDeltaLambdaBar(double d) { deltaLambdaBar = d; }
