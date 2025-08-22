// OriHinge.cpp
#include "OriHinge.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <Parameter.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include <ElementResponse.h>

#include <elementAPI.h>

// Nose
static int numOriHinge = 0;

#define OPS_Export
OPS_Export void *OPS_OriHinge()
{
	Element *theElement = 0;
	if (OPS_GetNumRemainingInputArgs() < 5)
	{
		opserr << "ERROR: insufficient args for OriHinge element\n";
		return 0;
	}

	int tag, nd1, nd2, nd3, nd4;
	int numArgs = 5;
	int iData[5];
	if (OPS_GetIntInput(&numArgs, iData) != 0)
		return nullptr;

	tag = iData[0];
	nd1 = iData[1];
	nd2 = iData[2];
	nd3 = iData[3];
	nd4 = iData[4];

	return new OriHinge(tag, nd1, nd2, nd3, nd4);
}

// ----- Constructores -----
OriHinge::OriHinge()
	: Element(0, ELE_TAG_OriHinge), connectedExternalNodes(4), theMatrix(0), theMass(0), theVector(0), J(12), d2thetadxi2(12, 12), theLoad(0)
{
	for (int i = 0; i < 4; i++)
		theNodes[i] = 0;
	theta0 = 0.0;
	theta = 0.0;
}

OriHinge::OriHinge(int tag, int node1, int node2, int node3, int node4)
	: Element(tag, ELE_TAG_OriHinge), connectedExternalNodes(4), theMatrix(0), theMass(0), theVector(0), J(12), d2thetadxi2(12, 12), theLoad(0)
{
	connectedExternalNodes(0) = node1;
	connectedExternalNodes(1) = node2;
	connectedExternalNodes(2) = node3;
	connectedExternalNodes(3) = node4;

	for (int i = 0; i < 4; i++)
		theNodes[i] = 0;
	theta0 = 0.0;
	theta = 0.0;
	Matrix K(24, 24);
	Matrix M(24, 24);
	Vector F(24);
	Vector J(12);
	Matrix d2thetadxi2(12, 12);
}

OriHinge::~OriHinge() {}

// ----- Métodos obligatorios -----
const ID &OriHinge::getExternalNodes()
{
	return connectedExternalNodes;
}

Node **OriHinge::getNodePtrs()
{
	return theNodes;
}

void OriHinge::setDomain(Domain *theDomain)
{
	if (theDomain == 0)
	{
		opserr << "OriHinge::setDomain - theDomain is null\n";
		return;
	}

	for (int i = 0; i < 4; i++)
	{
		theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
		if (theNodes[i] == 0)
		{
			opserr << "OriHinge::setDomain - node " << connectedExternalNodes(i) << " does not exist in the domain\n";
			return;
		}
	}

	for (int i = 0; i < 4; i++)
	{
		if (theNodes[i]->getNumberDOF() != 6)
		{
			opserr << "OriHinge::setDomain() - Node "
				   << connectedExternalNodes(i)
				   << " does not have 6 DOF.\n";
		}
	}

	this->DomainComponent::setDomain(theDomain);
	calculateVectors();
	theta0 = calculateTheta();
	theta = calculateTheta();
}

int OriHinge::commitState()
{
	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0)
	{
		opserr << "CorotTruss::commitState () - failed in base class\n";
	}
	// retVal = theMaterial->commitState();
	return retVal;
}

int OriHinge::revertToLastCommit()
{
	theta = theta0;
	return 0;
}

int OriHinge::revertToStart()
{
	theta0 = 0.0;
	theta = 0.0;
	return 0;
}

int OriHinge::update()
{
	calculateVectors();
	theta = calculateThetaFromU();
	return 0;
}

void OriHinge::calculateVectors()
{
	// TODO Calcular los vectores directores incluye desplazamientos
	// Calculates jacobian and d2thetadxi2
	J.Zero();
	d2thetadxi2.Zero();
}
float OriHinge::calculateTheta()
{
	// TODO Calcular el ángulo theta
	return 0.0;
}

float OriHinge::calculateThetaFromU()
{
	// TODO Calcular el ángulo theta a partir de los desplazamientos
	return 0.0;
}

float OriHinge::getMoment(float theta)
{
	// TODO Calcular el momento a partir del ángulo theta
	float k = 1e6; // Rigidez muy alta
	return k * (theta - theta0);
}

// Vector OriHinge::jacobian()
// {
// 	// Calcular la matriz jacobiana
// 	Vector J(12);
// 	J.Zero();
// 	return J;
// }

Matrix OriHinge::outer(Vector a, Vector b)
{
	Matrix m(a.Size(), b.Size());
	m.Zero();
	for (int i = 0; i < a.Size(); i++)
		for (int j = 0; j < b.Size(); j++)
			m(i, j) = a(i) * b(j);
	return m;
}

Matrix OriHinge::dia(Vector a, Vector b)
{
	return outer(a, b) + outer(b, a);
}

const Matrix &OriHinge::getTangentStiff()
{
	Matrix &K = *theMatrix;
	K.Zero();
	int ndof = 6;
	double k = 1e6; // Rigidez muy alta
	double moment = getMoment(theta);
	Matrix kg = moment * d2thetadxi2;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			K(i * ndof, j * ndof) = k * J(i * 3) * J(j * 3) + kg(i * 3, j * 3);
			K(i * ndof, j * ndof + 1) = k * J(i * 3) * J(j * 3 + 1) + kg(i * 3, j * 3 + 1);
			K(i * ndof, j * ndof + 2) = k * J(i * 3) * J(j * 3 + 2) + kg(i * 3, j * 3 + 2);

			K(i * ndof + 1, j * ndof) = k * J(i * 3 + 1) * J(j * 3) + kg(i * 3 + 1, j * 3);
			K(i * ndof + 1, j * ndof + 1) = k * J(i * 3 + 1) * J(j * 3 + 1) + kg(i * 3 + 1, j * 3 + 1);
			K(i * ndof + 1, j * ndof + 2) = k * J(i * 3 + 1) * J(j * 3 + 2) + kg(i * 3 + 1, j * 3 + 2);

			K(i * ndof + 2, j * ndof) = k * J(i * 3 + 2) * J(j * 3) + kg(i * 3 + 2, j * 3);
			K(i * ndof + 2, j * ndof + 1) = k * J(i * 3 + 2) * J(j * 3 + 1) + kg(i * 3 + 2, j * 3 + 1);
			K(i * ndof + 2, j * ndof + 2) = k * J(i * 3 + 2) * J(j * 3 + 2) + kg(i * 3 + 2, j * 3 + 2);
		}
	}

	return K;
}

const Matrix &OriHinge::getInitialStiff()
{
	Matrix &K = *theMatrix;
	K.Zero();
	int ndof = 6;
	double k = 1e6; // Rigidez muy alta
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			K(i * ndof, j * ndof) = J(i * 3) * J(j * 3);
			K(i * ndof, j * ndof + 1) = J(i * 3) * J(j * 3 + 1);
			K(i * ndof, j * ndof + 2) = J(i * 3) * J(j * 3 + 2);

			K(i * ndof + 1, j * ndof) = J(i * 3 + 1) * J(j * 3);
			K(i * ndof + 1, j * ndof + 1) = J(i * 3 + 1) * J(j * 3 + 1);
			K(i * ndof + 1, j * ndof + 2) = J(i * 3 + 1) * J(j * 3 + 2);

			K(i * ndof + 2, j * ndof) = J(i * 3 + 2) * J(j * 3);
			K(i * ndof + 2, j * ndof + 1) = J(i * 3 + 2) * J(j * 3 + 1);
			K(i * ndof + 2, j * ndof + 2) = J(i * 3 + 2) * J(j * 3 + 2);
		}
	}
	K *= k;
	return *theMatrix;
}

const Matrix &OriHinge::getMass()
{
	Matrix &KK = *theMass;
	KK.Zero();
	return *theMass;
}

const Vector &OriHinge::getResistingForce()
{
	Vector &F = *theVector;
	F.Zero();
	double moment = getMoment(theta);
	F = J * moment; // checks out
	return *theVector;
}

const Vector &OriHinge::getResistingForceIncInertia()
{
	Vector &F = *theVector;
	F.Zero();
	double moment = getMoment(theta);
	F = J * moment; // checks out
	return *theVector;
}

void OriHinge::Print(OPS_Stream &s, int flag)
{
	s << "OriHinge element, tag: " << this->getTag() << "\n";
	s << "Connected nodes: " << connectedExternalNodes;
}

const Matrix &OriHinge::getDamp(void)
{

	Matrix &K = *theMatrix;
	K.Zero();

	return *theMatrix;
}
void OriHinge::zeroLoad(void)
{
	theLoad->Zero();
}

int OriHinge::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	opserr << "CorotTruss::addLoad - load type unknown for OriHinge with tag: " << this->getTag() << endln;

	return -1;
}

int OriHinge::addInertiaLoadToUnbalance(const Vector &accel)
{

	return 0;
}

Response *
OriHinge::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "OriHinge");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);
	output.attr("node3", connectedExternalNodes[2]);
	output.attr("node4", connectedExternalNodes[3]);
	output.endTag();
	return theResponse;
}

int OriHinge::getResponse(int responseID, Information &eleInfo)
{
	double strain;

	return 0;
}
int OriHinge::setParameter(const char **argv, int argc, Parameter &param)
{
	return -1;
}

int OriHinge::updateParameter(int parameterID, Information &info)
{
	return -1;
}
int OriHinge::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int OriHinge::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

	return 0;
}