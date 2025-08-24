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
#include <cmath>
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
	double kf;
	numArgs = 1;
	if (OPS_GetDouble(&numArgs, &kf) != 0)
	{
		opserr << "WARNING invalid double kf in element OriHinge " << endln;
		return 0;
	}

	tag = iData[0];
	nd1 = iData[1];
	nd2 = iData[2];
	nd3 = iData[3];
	nd4 = iData[4];

	theElement = new OriHinge(tag, nd1, nd2, nd3, nd4, kf);
	if (theElement == 0)
	{
		opserr << "WARNING: out of memory: element OriHinge " << iData[0] << " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
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
	double kf = 0.0;

	static std::map<int, Vector> meshdata;
	if (info(0) == 1)
	{

		int numRemainingArgs = OPS_GetNumRemainingInputArgs();
		if (numRemainingArgs < 5)
		{
			opserr << "Invalid Args want: element OriHinge $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
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
		numData = 5;
		if (OPS_GetDouble(&numData, &kf) != 0)
		{
			opserr << "WARNING: Invalid kf: \n";
			return 0;
		}

		if (info.Size() < 2)
		{
			opserr << "WARNING: need info -- inmesh, meshtag\n";
			return 0;
		}

		// save the data for a mesh
		Vector &mdata = meshdata[info(1)];
		mdata.resize(5);
		mdata(0) = (double)nd1;
		mdata(1) = (double)nd1;
		mdata(2) = (double)nd3;
		mdata(3) = (double)nd4;
		mdata(4) = (double)kf;
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
		if (mdata.Size() < 5)
			return 0;

		iData[0] = info(2);
		iData[1] = info(3);
		iData[2] = info(4);
		iData[3] = info(5);
		iData[4] = info(6);
		kf = info(7);
	}
	theElement = new OriHinge(iData[0], iData[1], iData[2], iData[3], iData[4], kf);

	if (theElement == 0)
	{
		opserr << "WARNING: out of memory: element OriHinge " << iData[0] << " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
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
	kf = 0.0;
}

OriHinge::OriHinge(int tag, int node1, int node2, int node3, int node4, double pkf)
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
	kf = pkf;
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
		opserr << "OriHinge::commitState () - failed in base class\n";
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
	J.Zero();
	d2thetadxi2.Zero();

	Vector i = theNodes[0]->getCrds() + theNodes[0]->getTrialDisp();
	Vector j = theNodes[1]->getCrds() + theNodes[1]->getTrialDisp();
	Vector k = theNodes[2]->getCrds() + theNodes[2]->getTrialDisp();
	Vector l = theNodes[3]->getCrds() + theNodes[3]->getTrialDisp();

	// opserr << "Calculate vectors!" << endln;

	// Imprimir trial displacements

	// opserr << "Node 1 trial disp: " << theNodes[0]->getTrialDisp() << endln;
	// opserr << "Node 2 trial disp: " << theNodes[1]->getTrialDisp() << endln;
	// opserr << "Node 3 trial disp: " << theNodes[2]->getTrialDisp() << endln;
	// opserr << "Node 4 trial disp: " << theNodes[3]->getTrialDisp() << endln;

	Vector rij = i - j;
	Vector rkj = k - j;
	Vector rkl = k - l;

	// Normalizar el vector en nuevas variables rij_n
	double rij_n = rij.Norm();
	if (rij_n == 0)
	{
		opserr << "OriHinge::calculateVectors() - Zero length vector between nodes " << connectedExternalNodes(0) << " and " << connectedExternalNodes(1) << endln;
		return;
	}
	double rkj_n = rkj.Norm();
	if (rkj_n == 0)
	{
		opserr << "OriHinge::calculateVectors() - Zero length vector between nodes " << connectedExternalNodes(2) << " and " << connectedExternalNodes(1) << endln;
		return;
	}
	double rkl_n = rkl.Norm();
	if (rkl_n == 0)
	{
		opserr << "OriHinge::calculateVectors() - Zero length vector between nodes " << connectedExternalNodes(2) << " and " << connectedExternalNodes(3) << endln;
		return;
	}

	Vector m = cross(rij, rkj); // rij ^ rkj;
	double m_n2 = m.Norm() * m.Norm();
	Vector n = cross(rkj, rkl); // rkj ^ rkl;
	double n_n2 = n.Norm() * n.Norm();

	// A = rij@rkj/rkj_n**2
	double A = (rij ^ rkj) / (rkj_n * rkj_n);
	// B = rkl@rkj/self.rkj_n**2
	double B = (rkl ^ rkj) / (rkj_n * rkj_n);

	// opserr << todas las cantidades previamente calculadas
	// opserr << "rij: " << rij << endln;
	// opserr << "rkj: " << rkj << endln;
	// opserr << "rkl: " << rkl << endln;
	// opserr << "rij_n: " << rij_n << endln;
	// opserr << "rkj_n: " << rkj_n << endln;
	// opserr << "rkl_n: " << rkl_n << endln;
	// opserr << "m: " << m << endln;
	// opserr << "n: " << n << endln;
	// opserr << "A: " << A << endln;
	// opserr << "B: " << B << endln;

	// dA_dxj = (1 / rkj_n**2) *  ((2 * A - 1) * rkj - rij)
	Vector dA_dxj = (1.0 / (rkj_n * rkj_n)) * ((2.0 * A - 1.0) * rkj - rij);
	// dB_dxj = (1 / rkj_n**2) * (2 * B * rkj - rkl)
	Vector dB_dxj = (1.0 / (rkj_n * rkj_n)) * (2.0 * B * rkj - rkl);
	// self.dA_dxk = (1 / self.rkj_n**2) * (-2 * self.A * self.rkj + self.rij)
	Vector dA_dxk = (1.0 / (rkj_n * rkj_n)) * (-2.0 * A * rkj + rij);
	// self.dB_dxk = (1 / self.rkj_n**2) * ((1 - 2 * self.B) * self.rkj + self.rkl)
	Vector dB_dxk = (1.0 / (rkj_n * rkj_n)) * ((1.0 - 2.0 * B) * rkj + rkl);

	// opserr << "dA_dxj: " << dA_dxj << endln;
	// opserr << "dB_dxj: " << dB_dxj << endln;
	// opserr << "dA_dxk: " << dA_dxk << endln;
	// opserr << "dB_dxk: " << dB_dxk << endln;

	// self.dtdxi = (self.rkj_n/self.m_n2) * self.m
	Vector dtdxi = (rkj_n / m_n2) * m;
	// self.dtdxl = (-self.rkj_n/self.n_n2) * self.n
	Vector dtdxl = (-rkj_n / n_n2) * n;
	// self.dtdxj = (np.dot(self.rij, self.rkj)/self.rkj_n**2-1)*self.dtdxi - np.dot(self.rkl, self.rkj)/self.rkj_n**2 * self.dtdxl
	Vector dtdxj = (A - 1.0) * dtdxi - B * dtdxl;
	// self.dtdxk = (np.dot(self.rkl, self.rkj)/self.rkj_n**2-1)*self.dtdxl - np.dot(self.rij, self.rkj)/self.rkj_n**2 * self.dtdxi
	Vector dtdxk = (B - 1.0) * dtdxl - A * dtdxi;

	// opserr << "dtdxi: " << dtdxi << endln;
	// opserr << "dtdxl: " << dtdxl << endln;
	// opserr << "dtdxj: " << dtdxj << endln;
	// opserr << "dtdxk: " << dtdxk << endln;

	// cross_product = np.cross(self.rkj, self.m)
	//  d2tdxi2 = -(self.rkj_n / (self.m_n2**2)) * self.dia(self.m, cross_product)  # yes
	Vector cross_product = cross(rkj, m);
	Matrix d2tdxi2 = -(rkj_n / (m_n2 * m_n2)) * dia(m, cross_product);

	// opserr << "d2tdxi2: " << d2tdxi2 << endln;

	// cross_product_n = np.cross(self.rkj, self.n)
	// d2tdxl2 = (self.rkj_n / (self.n_n2**2)) * self.dia(self.n, cross_product_n)  # yes
	cross_product = cross(rkj, n);
	Matrix d2tdxl2 = (rkj_n / (n_n2 * n_n2)) * dia(n, cross_product);

	// opserr << "d2tdxl2: " << d2tdxl2 << endln;

	// cross_product_m_rij = np.cross(self.rij, self.m)
	// d2tdxixk = (np.outer(self.m, self.rkj) / (self.m_n2 * self.rkj_n)) + (self.rkj_n / (self.m_n2**2)) * self.dia(self.m, cross_product_m_rij)  # yes
	cross_product = cross(rij, m);
	Matrix d2tdxixk = (outer(m, rkj) / (m_n2 * rkj_n)) + (rkj_n / (m_n2 * m_n2)) * dia(m, cross_product);

	// opserr << "d2tdxixk: " << d2tdxixk << endln;

	// cross_product_n_rkl = np.cross(self.rkl, self.n)
	// d2tdxlxj = (np.outer(self.n, self.rkj) / (self.n_n2 * self.rkj_n)) - (self.rkj_n / (self.n_n2**2)) * self.dia(self.n, cross_product_n_rkl)  # yes
	cross_product = cross(rkl, n);
	Matrix d2tdxlxj = (outer(n, rkj) / (n_n2 * rkj_n)) - (rkj_n / (n_n2 * n_n2)) * dia(n, cross_product);

	// opserr << "d2tdxlxj: " << d2tdxlxj << endln;

	// cross_product_m_diff = np.cross(self.rkj - self.rij, self.m)
	// d2tdxixj = -(np.outer(self.m, self.rkj) / (self.m_n2 * self.rkj_n)) + (self.rkj_n / (self.m_n2**2)) * self.dia(self.m, cross_product_m_diff)  # yes
	cross_product = cross(rkj - rij, m);
	Matrix d2tdxixj = -1.0 * (outer(m, rkj) / (m_n2 * rkj_n)) + (rkj_n / (m_n2 * m_n2)) * dia(m, cross_product);

	// opserr << "d2tdxixj: " << d2tdxixj << endln;

	// cross_product_n_diff = np.cross(self.rkj - self.rkl, self.n)
	// d2tdxlxk = -(np.outer(self.n, self.rkj) / (self.n_n2 * self.rkj_n)) - (self.rkj_n / (self.n_n2**2)) * self.dia(self.n, cross_product_n_diff)  # yes
	cross_product = cross(rkj - rkl, n);
	Matrix d2tdxlxk = -1.0 * (outer(n, rkj) / (n_n2 * rkj_n)) - (rkj_n / (n_n2 * n_n2)) * dia(n, cross_product);

	// opserr << "d2tdxlxk: " << d2tdxlxk << endln;

	// d2tdxj2 = np.outer(self.dtdxi, self.dA_dxj) + (self.A - 1) * d2tdxixj - (np.outer(self.dtdxl, self.dB_dxj) + self.B * d2tdxlxj)
	Matrix d2tdxj2 = outer(dtdxi, dA_dxj) + (A - 1.0) * d2tdxixj - (outer(dtdxl, dB_dxj) + B * d2tdxlxj);

	// opserr << "d2tdxj2: " << d2tdxj2 << endln;

	// d2tdxjxk = np.outer(self.dtdxi, self.dA_dxk) + (self.A - 1) * d2tdxixk - (np.outer(self.dtdxl, self.dB_dxk) + self.B * d2tdxlxk)
	Matrix d2tdxjxk = outer(dtdxi, dA_dxk) + (A - 1.0) * d2tdxixk - (outer(dtdxl, dB_dxk) + B * d2tdxlxk);

	// opserr << "d2tdxjxk: " << d2tdxjxk << endln;

	// d2tdxk2 = np.outer(self.dtdxl, self.dB_dxk) + (self.B - 1) * d2tdxlxk - (np.outer(self.dtdxi, self.dA_dxk) + self.A * d2tdxixk)
	Matrix d2tdxk2 = outer(dtdxl, dB_dxk) + (B - 1.0) * d2tdxlxk - (outer(dtdxi, dA_dxk) + A * d2tdxixk);

	// opserr << "d2tdxk2: " << d2tdxk2 << endln;

	// d2tdxlxi = np.zeros((3, 3))
	Matrix d2tdxlxi(3, 3);
	d2tdxlxi.Zero();
	// d2theta_dxi2 = np.block([[d2tdxi2, d2tdxixj, d2tdxixk, d2tdxlxi],
	//                                   [d2tdxixj.T, d2tdxj2, d2tdxjxk, d2tdxlxj.T],
	//                                   [d2tdxixk.T, d2tdxjxk.T,
	//                                       d2tdxk2, d2tdxlxk.T],
	//                                   [d2tdxlxi.T, d2tdxlxj, d2tdxlxk, d2tdxl2]])
	for (int r = 0; r < 3; r++)
	{

		for (int c = 0; c < 3; c++)
		{
			d2thetadxi2(r, c) = d2tdxi2(r, c);
			d2thetadxi2(r, c + 3) = d2tdxixj(r, c);
			d2thetadxi2(r, c + 6) = d2tdxixk(r, c);
			d2thetadxi2(r, c + 9) = d2tdxlxi(r, c);

			d2thetadxi2(r + 3, c) = d2tdxixj(c, r);
			d2thetadxi2(r + 3, c + 3) = d2tdxj2(r, c);
			d2thetadxi2(r + 3, c + 6) = d2tdxjxk(r, c);
			d2thetadxi2(r + 3, c + 9) = d2tdxlxj(c, r);

			d2thetadxi2(r + 6, c) = d2tdxixk(c, r);
			d2thetadxi2(r + 6, c + 3) = d2tdxjxk(c, r);
			d2thetadxi2(r + 6, c + 6) = d2tdxk2(r, c);
			d2thetadxi2(r + 6, c + 9) = d2tdxlxk(c, r);

			d2thetadxi2(r + 9, c) = d2tdxlxi(c, r);
			d2thetadxi2(r + 9, c + 3) = d2tdxlxj(r, c);
			d2thetadxi2(r + 9, c + 6) = d2tdxlxk(r, c);
			d2thetadxi2(r + 9, c + 9) = d2tdxl2(r, c);
		}
	}

	// opserr << "d2thetadxi2: " << d2thetadxi2 << endln;

	// J = np.array( [self.dtdxi, self.dtdxj, self.dtdxk, self.dtdxl])

	J(0) = dtdxi(0);
	J(1) = dtdxi(1);
	J(2) = dtdxi(2);
	J(3) = dtdxj(0);
	J(4) = dtdxj(1);
	J(5) = dtdxj(2);
	J(6) = dtdxk(0);
	J(7) = dtdxk(1);
	J(8) = dtdxk(2);
	J(9) = dtdxl(0);
	J(10) = dtdxl(1);
	J(11) = dtdxl(2);

	// opserr << "J: " << J << endln;
}

Vector OriHinge::cross(const Vector &a, const Vector &b)
{
	if (a.Size() != 3 || b.Size() != 3)
	{
		opserr << "OriHinge::cross() - Vectors must be of size 3\n";
		Vector c(3);
		c.Zero();
		return c;
	}
	Vector c(3);
	c(0) = a(1) * b(2) - a(2) * b(1);
	c(1) = a(2) * b(0) - a(0) * b(2);
	c(2) = a(0) * b(1) - a(1) * b(0);
	return c;
}

double OriHinge::calculateTheta()
{
	// eta = 1
	//  if round(np.dot(self.m, self.rkl), 8) != 0:
	//      eta = np.sign(np.dot(self.m, self.rkl))
	//  cosa  = np.dot(self.m, self.n) / \
	//     (np.linalg.norm(self.m)*np.linalg.norm(self.n))
	//  cosa = max(min(cosa, 1), -1)
	//  theta_ = eta*np.arccos(cosa)
	//  theta_ = theta_ % (2*np.pi)

	// return theta_
	// Ahora en C++

	double eta = 1.0;
	Vector i = theNodes[0]->getCrds() + theNodes[0]->getTrialDisp();
	Vector j = theNodes[1]->getCrds() + theNodes[1]->getTrialDisp();
	Vector k = theNodes[2]->getCrds() + theNodes[2]->getTrialDisp();
	Vector l = theNodes[3]->getCrds() + theNodes[3]->getTrialDisp();
	Vector rij = i - j;
	Vector rkj = k - j;
	Vector rkl = k - l;
	Vector m = cross(rij, rkj); // rij ^ rkj;
	Vector n = cross(rkj, rkl); // rkj ^ rkl;
	double m_n = m.Norm();
	double n_n = n.Norm();
	double dot_m_rkl = m ^ rkl;

	if (!(abs(dot_m_rkl) < 1e-6))
	{
		eta = (dot_m_rkl > 0) ? 1.0 : -1.0;
	}
	double cosa = (m ^ n) / (m_n * n_n);
	if (cosa > 1.0)
		cosa = 1.0;
	if (cosa < -1.0)
		cosa = -1.0;
	double theta_ = eta * acos(cosa);
	theta_ = fmod(theta_, (2.0 * PI)); // Esta función es distinta que la de Python
	if (theta_ < 0)
	{
		theta_ += 2.0 * PI;
	}
	return theta_;
}

double OriHinge::calculateThetaFromU()
{

	// Ue = self.Ue.T.flatten().reshape([12, 1])
	// J = self.jacobian().flatten().reshape([12, 1])
	// r = self.theta_0 + J.T@Ue
	// r = r % (2*np.pi)

	// return r
	Vector Ue(12);
	for (int i = 0; i < 4; i++)
	{
		Vector disp = theNodes[i]->getTrialDisp();
		Ue(i * 3) = disp(0);
		Ue(i * 3 + 1) = disp(1);
		Ue(i * 3 + 2) = disp(2);
	}
	double r = theta0 + (J ^ Ue);
	r = fmod(r, (2.0 * PI));
	if (r < 0)
	{
		r += 2.0 * PI;
	}
	return r;
}

double OriHinge::getMoment(double ptheta)
{
	double Lr = 0.0;
	// Calcular Lr como norma de rkj
	Vector j = theNodes[1]->getCrds() + theNodes[1]->getTrialDisp();
	Vector k = theNodes[2]->getCrds() + theNodes[2]->getTrialDisp();
	Vector rkj = k - j;
	Lr = rkj.Norm();
	return Lr * kf * (ptheta - theta0);
}

double OriHinge::getKf(double ptheta)
{
	double Lr = 0.0;
	Vector j = theNodes[1]->getCrds() + theNodes[1]->getTrialDisp();
	Vector k = theNodes[2]->getCrds() + theNodes[2]->getTrialDisp();
	Vector rkj = k - j;
	Lr = rkj.Norm();
	return Lr * kf;
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
	opserr << "Entra a la matrix tangete\n";
	opserr << "Theta " << theta * 180.0 / PI << endln;
	opserr << "Theta_0 " << theta0 * 180.0 / PI << endln;
	Matrix K = theMatrix;
	int ndof = 6;
	double k = getKf(theta); // Rigidez muy alta
	double moment = getMoment(theta);
	opserr << "M " << moment << endln;
	Matrix kg = moment * d2thetadxi2;
	Matrix ke = outer(J, J) * k;
	Matrix kt = ke + kg;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			K(i * ndof, j * ndof) = kt(i * 3, j * 3);
			K(i * ndof, j * ndof + 1) = kt(i * 3, j * 3 + 1);
			K(i * ndof, j * ndof + 2) = kt(i * 3, j * 3 + 2);

			K(i * ndof + 1, j * ndof) = kt(i * 3 + 1, j * 3);
			K(i * ndof + 1, j * ndof + 1) = kt(i * 3 + 1, j * 3 + 1);
			K(i * ndof + 1, j * ndof + 2) = kt(i * 3 + 1, j * 3 + 2);

			K(i * ndof + 2, j * ndof) = kt(i * 3 + 2, j * 3);
			K(i * ndof + 2, j * ndof + 1) = kt(i * 3 + 2, j * 3 + 1);
			K(i * ndof + 2, j * ndof + 2) = kt(i * 3 + 2, j * 3 + 2);
		}
	}

	// opserr << "K: " << K << endln;

	return theMatrix;
}

const Matrix &OriHinge::getInitialStiff()
{
	opserr << "Entra a la matrix inicial\n";
	opserr << "Theta " << theta * 180.0 / PI << endln;
	opserr << "Theta_0 " << theta0 * 180.0 / PI << endln;
	Matrix K = theMatrix;
	int ndof = 6;
	double k = getKf(theta); // Rigidez muy alta
	double moment = getMoment(theta);
	opserr << "M " << moment << endln;
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
	opserr << F << endln;
	opserr << moment << endln;
	opserr << J << endln;

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
	opserr << "OriHinge::addLoad - load type unknown for OriHinge with tag: " << this->getTag() << endln;

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