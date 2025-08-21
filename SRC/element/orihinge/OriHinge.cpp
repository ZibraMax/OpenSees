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

Matrix OriHinge::K(12, 12);
Matrix OriHinge::M(12, 12);
Vector OriHinge::F(12);

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
	: Element(0, ELE_TAG_OriHinge), connectedExternalNodes(4)
{
	for (int i = 0; i < 4; i++)
		theNodes[i] = 0;
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
		if (theNodes[i]->getNumberDOF() != 3)
		{
			opserr << "OriHinge::setDomain() - Node "
				   << connectedExternalNodes(i)
				   << " does not have 3 DOF.\n";
		}
	}

	this->DomainComponent::setDomain(theDomain);
}

const Matrix &OriHinge::getTangentStiff()
{
	K.Zero();
	return K;
}

const Matrix &OriHinge::getInitialStiff()
{
	K.Zero();
	return K;
}

const Matrix &OriHinge::getMass()
{
	M.Zero();
	return M;
}

const Vector &OriHinge::getResistingForce()
{
	F.Zero();
	return F;
}

const Vector &OriHinge::getResistingForceIncInertia()
{
	F.Zero();
	return F;
}

void OriHinge::Print(OPS_Stream &s, int flag)
{
	s << "OriHinge element, tag: " << this->getTag() << "\n";
	s << "Connected nodes: " << connectedExternalNodes;
}
