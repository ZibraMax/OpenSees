// OriHinge.cpp
#include "OriHinge.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>

#include <Parameter.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include <ElementResponse.h>

Matrix OriHinge::theMatrix(24, 24);
Matrix OriHinge::theMass(24, 24);
Vector OriHinge::theVector(12);
Vector OriHinge::theLoad(12);
Vector OriHinge::J(12);
Matrix OriHinge::d2thetadxi2(12, 12);
#include <elementAPI.h>

#define OPS_Export
OPS_Export void *OPS_OriHingeElement()
{
	Element *theElement = 0;
	int numRemainingArgs = OPS_GetNumRemainingInputArgs();
	if (numRemainingArgs < 5)
	{
		opserr << "ERROR: insufficient args for OriHinge element\n";
		return 0;
	}

	int tag, nd1, nd2, nd3, nd4;
	int iData[5];
	int numArgs = 5;
	if (OPS_GetInt(&numArgs, iData) != 0)
	{

		opserr << "WARNING invalid integer (tag, nd1, nd2, nd3, nd4) in element OriHinge " << endln;
		return 0;
	}

	tag = iData[0];
	nd1 = iData[1];
	nd2 = iData[2];
	nd3 = iData[3];
	nd4 = iData[4];
	theElement = new OriHinge(tag, nd1, nd2, nd3, nd4);
	if (theElement == 0)
	{
		opserr << "WARNING: out of memory: element CorotTruss " << iData[0] << " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	}

	return theElement;
}

void *
OPS_OriHingeElement(const ID &info)
{
	if (info.Size() == 0)
		return 0;

	Element *theElement = 0;

	int iData[5];

	int nd1 = 0;
	int nd2 = 0;
	int nd3 = 0;
	int nd4 = 0;

	static std::map<int, Vector> meshdata;
	if (info(0) == 1)
	{

		int numRemainingArgs = OPS_GetNumRemainingInputArgs();
		if (numRemainingArgs < 4)
		{
			opserr << "Invalid Args want: element CorotTruss $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
			return 0;
		}
		int numData = 1;
		if (OPS_GetInt(&numData, &nd1) != 0)
		{
			opserr << "WARNING: Invalid nd1: \n";
			return 0;
		}
		numData = 2;
		if (OPS_GetInt(&numData, &nd2) != 0)
		{
			opserr << "WARNING: Invalid nd2: \n";
			return 0;
		}
		numData = 3;
		if (OPS_GetInt(&numData, &nd3) != 0)
		{
			opserr << "WARNING: Invalid nd3: \n";
			return 0;
		}
		numData = 4;
		if (OPS_GetInt(&numData, &nd4) != 0)
		{
			opserr << "WARNING: Invalid nd4: \n";
			return 0;
		}

		if (info.Size() < 2)
		{
			opserr << "WARNING: need info -- inmesh, meshtag\n";
			return 0;
		}

		// save the data for a mesh
		Vector &mdata = meshdata[info(1)];
		mdata.resize(4);
		mdata(0) = (double)nd1;
		mdata(1) = (double)nd1;
		mdata(2) = (double)nd3;
		mdata(3) = (double)nd4;
		return &meshdata;
	}
	else if (info(0) == 2)
	{

		if (info.Size() < 7)
		{
			opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2\n";
			return 0;
		}

		Vector &mdata = meshdata[info(1)];
		if (mdata.Size() < 4)
			return 0;

		iData[0] = info(2);
		iData[1] = info(3);
		iData[2] = info(4);
		iData[3] = info(5);
		iData[4] = info(6);
	}
	theElement = new OriHinge(iData[0], iData[1], iData[2], iData[3], iData[4]);

	if (theElement == 0)
	{
		opserr << "WARNING: out of memory: element CorotTruss " << iData[0] << " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	}

	return theElement;
}

// ----- Constructores -----
OriHinge::OriHinge()
	: Element(0, ELE_TAG_OriHinge), connectedExternalNodes(4)
{
	for (int i = 0; i < 4; i++)
		theNodes[i] = 0;
	theta0 = 0.0;
	theta = 0.0;
}

OriHinge::OriHinge(int tag, int node1, int node2, int node3, int node4)
	: Element(tag, ELE_TAG_OriHinge), connectedExternalNodes(4)
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
	opserr << "Entra a domain\n";
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
	J += 1;
	d2thetadxi2.Zero();
	d2thetadxi2 += 1;
}
float OriHinge::calculateTheta()
{
	// TODO Calcular el ángulo theta
	return 30.0 * 180.0 / PI;
}

float OriHinge::calculateThetaFromU()
{
	// TODO Calcular el ángulo theta a partir de los desplazamientos
	return 32.0 * 180.0 / PI;
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
	opserr << "EntrL a la matrix\n";
	Matrix K = theMatrix;
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

	return theMatrix;
}

const Matrix &OriHinge::getInitialStiff()
{
	opserr << "EntrL a la matrix init\n";
	Matrix K = theMatrix;
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
	return theMatrix;
}

const Matrix &OriHinge::getMass()
{
	Matrix KK = theMass;
	KK.Zero();
	return theMass;
}

const Vector &OriHinge::getResistingForce()
{
	Vector F = theVector;
	F.Zero();
	double moment = getMoment(theta);
	F = J * moment; // checks out
	return theVector;
}

const Vector &OriHinge::getResistingForceIncInertia()
{
	Vector F = theVector;
	F.Zero();
	double moment = getMoment(theta);
	F = J * moment; // checks out
	return theVector;
}

void OriHinge::Print(OPS_Stream &s, int flag)
{
	s << "OriHinge element, tag: " << this->getTag() << "\n";
	s << "Connected nodes: " << connectedExternalNodes;
}

const Matrix &OriHinge::getDamp(void)
{

	Matrix K = theMatrix;
	K.Zero();

	return theMatrix;
}
void OriHinge::zeroLoad(void)
{
	theLoad.Zero();
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