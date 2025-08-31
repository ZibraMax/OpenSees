#ifndef MGDCM_H
#define MGDCM_H

#include <StaticIntegrator.h>
#include <Vector.h>

class LinearSOE;
class AnalysisModel;

class MGDCM : public StaticIntegrator
{
public:
	MGDCM(double deltaLambdaBar, int maxIter, int minIter, bool momentum);

	// core OpenSees interface
	int newStep(void);
	int update(const Vector &deltaU);
	int domainChanged(void);
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	// additional utilities
	void setDeltaLambdaBar(double d);

protected:
	double deltaLambdaBar;
	int increments;
	int maxIterMomentum;
	int minIterMomentum;
	bool useMomentum;

	double currentLambda;
	double dLambda;
	double signLastStep;

	Vector *duHat; // displacement due to reference load
	Vector *duBar; // correction vector
	Vector *du;	   // total increment
	Vector *phat;  // reference load vector
};

#endif
