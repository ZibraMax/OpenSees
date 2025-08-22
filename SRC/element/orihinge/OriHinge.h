#ifndef OriHinge_h
#define OriHinge_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

const double PI = 3.14159265358979323846;

class Node;
class Channel;

class OriHinge : public Element
{
public:
  OriHinge(int tag, int node1,
           int node2, int node3, int node4);

  OriHinge();
  ~OriHinge();

  const char *getClassType(void) const { return "OriHinge"; };

  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes(void) const { return 4; }
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);

  int getNumDOF(void) { return 24; } // 4 nodos × 6 GDL
  void setDomain(Domain *theDomain);

  // public methods to set the state of the element
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);
  void calculateVectors(void);
  float calculateTheta(void);
  float calculateThetaFromU(void);
  float getMoment(float theta);
  // Vector jacobian(void);
  Matrix dia(Vector a, Vector b);
  Matrix outer(Vector a, Vector b);

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getDamp(void);
  const Matrix &getMass(void);

  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  // int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);
  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

protected:
private:
  ID connectedExternalNodes; // contains the tags of the end nodes
  Node *theNodes[4];
  double theta0;
  double theta;
  double theta1 = 0.0;
  double theta2 = 0.0;

  static Vector theLoad;
  static Matrix theMass;
  static Matrix theMatrix;
  static Vector theVector;
  static Vector J;
  static Matrix d2thetadxi2;
};

#endif
