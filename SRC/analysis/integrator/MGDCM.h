#ifndef MGDCM_H
#define MGDCM_H

#include <StaticIntegrator.h>
#include <Vector.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Domain;

class MGDCM : public StaticIntegrator
{
public:
	MGDCM(double deltaLambdaBar, int maxIter, int minIter, bool momentum);
	~MGDCM();

	int newStep(void);
	int update(const Vector &deltaU);
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag = 0);

private:
	void getFext();
	double getDeltaLambda(Vector dup, Vector dur);
	double deltaLambdaBar;
	double numgsp;
	double lambda;

	int maxIterMomentum;
	int minIterMomentum;

	int i;
	int k;
	int sign;

	bool useMomentum;
	double dLambda;

	Vector *dupp1;
	Vector *dupc1;

	Vector *duHat;
	Vector *duBar;
	Vector *du;
	Vector *Fext;
};

#endif
