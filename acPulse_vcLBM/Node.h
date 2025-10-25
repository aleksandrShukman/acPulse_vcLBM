#pragma once

#include "Control.h"
#include "Vector3D.H"
#include <vector>
#include "Matrix.h"

using namespace std;

class Node
{
public:
	//!constructors
	Node() {} //!standard

	Node(double x, double y, double z, int indexI, int indexJ, int indexK, int level) : x_(x), y_(y), z_(z), indexI_(indexI), indexJ_(indexJ), indexK_(indexK), level_(level) //!init all data
	{
		for (long long dir = 0; dir < 19; dir++)
		{
			distr_[false][dir] = 0;
			distr_[true][dir]  = 0;

			distrPreCol_[dir] = 0;
		}
		turbVisc_ = 0;
	}

	//!destructor
	virtual ~Node()
	{

	}

	string tag_;
	int ID_;

	//!get and set methods
	virtual inline int    getIndexI() const { return indexI_; }
	virtual inline int    getIndexJ() const { return indexJ_; }
	virtual inline int    getIndexK() const { return indexK_; }
	virtual inline int    getLevel() const { return level_; }
    virtual inline int    getNormalDir() { return normalDir_; }
    virtual inline int    getNodeID() { return ID_; }
    virtual inline bool   getIsCorner() { return isCorner_; }
	virtual inline double getX() const { return x_; }
	virtual inline double getY() const { return y_; }
	virtual inline double getZ() const { return z_; }
	virtual inline double getNodeDensity() const { return densityNode_; }
	virtual inline double getNodeVelocityComponent(int dir) const { return velocityNode_[dir]; }
	virtual inline double getCollisionFrequency() const { return colFreq_; }
	virtual inline double getTurbulentViscosity() const { return turbVisc_; }
	virtual inline double getDistribution(bool& distr, int& dir) { if (dir < 19) return distr_[distr][dir]; else return NULL; };
	virtual inline double getCubicMachCorrection(int& dir) { if (dir < 19) return psi_[dir]; else return NULL; };
	virtual inline double getCubicMachCorrectionPreCol(int& dir) { if (dir < 19) return psiPreCol_[dir]; else return NULL; };
	virtual inline double getCubicMachCorrectionHalf(int& dir) { if (dir < 19) return psiHalf_[dir]; else return NULL; };
	virtual inline double getCubicMachCorrectionFull(int& dir) { if (dir < 19) return psiFull_[dir]; else return NULL; };
	virtual inline double getDistributionPreCol(int& dir) { if (dir < 19) return distrPreCol_[dir]; else return NULL; };
	virtual inline double getDistributionFictitious(int& dir) { if (dir < 19) return distrFictitious_[dir]; else return NULL; };
	virtual inline double getEquilibriumDistribution(int& dir, double* equil) { return equil[dir]; };
	virtual inline double getRROffEquilibrium(int& dir, double* nEquilRR) { return nEquilRR[dir]; };
	virtual inline double getStrainRateTensorFD(int i, int j) { return strainRateTensorFD_[i][j]; }
	virtual inline double getStrainRateTensorFDfull(int i, int j) { return strainRateTensorFDfull_[i][j]; }
	virtual inline double getStrainRateTensorFDhalf(int i, int j) { return strainRateTensorFDhalf_[i][j]; }
	virtual inline double getPi1Tensor(int i, int j) { return Pi1Tensor_[i][j]; }
	virtual inline double getTimeStep() const { return timeStep_; }
	virtual inline string getGridCouplingTag() { cerr << "\nGetGridCouplinTag() from node invoked." << flush; exit(1); } //!salvation implementation
    virtual inline string getTag() { return tag_; }
    virtual inline Node*  getNeighbor(int dir) { if (dir < 18) return neighbors_[dir]; else return NULL; }
    virtual inline Node*  getNormalNeighbor() { return normalNeighbor_; }
    virtual inline Node*  getTangentNeighbor1(int dir) { return tangentNeighbor1_[dir]; }
    virtual inline Node*  getTangentNeighbor2(int dir) { return tangentNeighbor2_[dir]; }
    virtual inline Node*  getNormalNeighborSnd() { return normalNeighborSnd_; }
    virtual inline Node*  getPartner() { return partner_; }

	virtual inline void	setDistributionPreCol(int& dir, double& distrVal) { if (dir < 19) distrPreCol_[dir] = distrVal; }
	virtual inline void	setCubicMachCorrection(int& dir, double psiVal) { if (dir < 19) psi_[dir] = psiVal; }
	virtual inline void	setCubicMachCorrectionPreCol(int& dir, double psiVal) { if (dir < 19) psiPreCol_[dir] = psiVal; }
	virtual inline void	setCubicMachCorrectionHalf(int& dir, double psiVal) { if (dir < 19) psiHalf_[dir] = psiVal; }
	virtual inline void	setCubicMachCorrectionFull(int& dir, double psiVal) { if (dir < 19) psiFull_[dir] = psiVal; }
	virtual inline void	setDistributionFictitious(int& dir, double& distrVal) { if (dir < 19) distrFictitious_[dir] = distrVal; }
	virtual inline void setStrainRateTensorFD(int i, int j, double value) { strainRateTensorFD_[i][j] = value; }
	virtual inline void setStrainRateTensorFDfull(int i, int j, double value) { strainRateTensorFDfull_[i][j] = value; }
	virtual inline void setStrainRateTensorFDhalf(int i, int j, double value) { strainRateTensorFDhalf_[i][j] = value; }
	virtual inline void setPi1Tensor(int i, int j, double value) { Pi1Tensor_[i][j] = value; }
	virtual inline void setDistribution(bool alternating, int& dir, double distrVal) { if (dir < 19) distr_[alternating][dir] = distrVal; }
	virtual inline void setCollisionFrequency(double& colFreq) { colFreq_ = colFreq;}
	virtual inline void setTurbulentViscosity(double& turbVisc) { turbVisc_ = turbVisc; }
	virtual inline void setLevel(int& level) { level_ = level; }
	virtual inline void setNeighbor(int& dir, Node*& neighbor) {if (dir<18) neighbors_[dir] = neighbor;}
	virtual inline void setPartner(Node*& partner) { partner_ = partner; }
	virtual inline void setNormalNeighbor(Node*& normalNeighbor) { normalNeighbor_ = normalNeighbor; }
	virtual inline void setNormalNeighborSnd(Node*& normalNeighborSnd) { normalNeighborSnd_ = normalNeighborSnd; }
	virtual inline void setNormalDir(int& dir) { normalDir_ = dir; }
	virtual inline void setTangentNeighbor1(Node*& tangentNeighbor1, int dir) { tangentNeighbor1_[dir] = tangentNeighbor1; }
	virtual inline void setTangentNeighbor2(Node*& tangentNeighbor2, int dir) { tangentNeighbor2_[dir] = tangentNeighbor2; }
	virtual inline void setIsCorner(bool value) { isCorner_ = value; }
	virtual inline void setTimeStep(double& timeStep) { timeStep_ = timeStep; }
	virtual inline void setNodeDensity(double& density) { densityNode_ = density; }
	virtual inline void setNodeVelocityComponent(int& dir, double& velocity) { velocityNode_[dir] = velocity; }

	//!members
	virtual void collideSRT(Control& bc, bool& alternating);
	virtual void collideRR(Control& bc, bool& alternating);
	virtual void collideMRT_GSbasis(Control& bc, bool& alternating);
	virtual void collideMRT_RMbasis(Control& bc, bool& alternating);
	virtual void collideHRR(Control& bc, bool& alternating);
	virtual void calcVelocity(Control& bc, bool& distr, double& dens, double* velo);
	virtual void calcVelocityPreCol(Control& bc, double& dens, double* velo);
	virtual void calcVelocityFictitious(Control& bc, double& dens, double* velo);
    virtual void calcEquilibrium2_GHbasis(double* equil, Control& bc, double& dens, double* velo);
    virtual void calcEquilibrium3_GHbasis(double* equil, Control& bc, double& dens, double* velo);
    virtual void calcEquilibrium4_GHbasis(double* equil, Control& bc, double& dens, double* velo);
    virtual void calcRRoffEquilibrium2_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo);
    virtual void calcRRoffEquilibrium3_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo);
    virtual void calcRRoffEquilibrium4_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo);
    virtual void calcHRRoffEquilibrium2_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo);
    virtual void calcHRRoffEquilibrium3_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo);
    virtual void calcHRRoffEquilibrium4_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo);
	virtual void calcPi1Tensor(bool& alternating, double* equil, Control& bc, double* velo);
	virtual void transformMomentsToDistributions_RMbasis(double* distributions, double* moments, Control& bC);
	virtual void calcAllMoments_GSbasis(Control& bc, double* moments, bool& alternating);
	virtual void calcAllMoments_GSbasis(Control& bc, double* distr, double* moments, bool& alternating);
	virtual void calcAllMoments_RMbasis(Control& bc, double* moments, double* distribution) const;
	virtual void calcEquilibrium4_RMbasis(double* distributions, double density, double* velo, Control& bc);
	virtual void calcEquilibrium4_RMbasis(double* distributions, double density, Vector3D* velo, Control& bc);
	virtual void calcAllMomentsPreCol_GSbasis(Control& bc, double* moments);
	virtual void calcEquilibriumMoments2_GSbasis(double* equilibriumMoments, double& density, double& momentumX, double& momentumY, double& momentumZ);
	virtual void transformMomentsToDistributions_GSbasis(double* distributions, double* moments, Control& bc);
	virtual void calcStrainRateTensorFDandCubicMachCorrection(Control& bc);
	virtual void calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(bool& alternatingL0, bool& alternatingL1, Control& bc);
	virtual void calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(bool& alternatingL0, bool& alternatingL1, Control& bc);
	virtual void calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(bool& alternatingL0, bool& alternatingL1, Control& bc);
	virtual void calcStrainRateTensorFDandCubicMachCorrectionFinePartners(bool& alternatingL1, Control& bc);
	virtual double calcDensity(Control& bc, bool& alternating);
	virtual double calcDensityPreCol(Control& bc);
	virtual double calcDensityFictitious();
	virtual double calcPressure(bool& alternating, Control& bc);
	virtual inline void transport(bool& /*distr*/) { cerr << "\nTransport from node invoked." << flush; exit(1); } //!salvation implementation
	virtual inline void transportPull(bool& /*distr*/) { cerr << "\nTransport (pull) from node invoked." << flush; exit(1); } //!salvation implementation

	//!inlet and outlet boundary conditions
	virtual void pressureOutlet(Control& bc, bool& alternating);
	virtual void LODIpressureOutlet(Control& bc, bool& alternating);
	virtual void LODIpressureCalcVelocity(Control& bc, double& domainLength, double& density, double& densityBoundaryNode, double& densityGradient, double* velocityGradient, double* veloBoundaryNode, double* velocity, bool negNormal, double* densityGradientT, double velocityGradientT[][3]);

	//!bounce-back methods
	virtual inline void halfWayBounceBack(bool &/*alternating*/, Control &/*bc*/) {}; //!virtual declaration for wallnode
	virtual inline void frictionless1stOrder(bool& /*alternating*/, Control& /*bc*/) {}; //!virtual declaration for wallnode

	//!new members from GridCoupling
	virtual inline void directCouplingCalcDeltaX(double* /*X*/, double* /*distr*/, double& /*rho*/, double* /*velo*/, string& /*interfaceType*/, Control& /*bc*/) { cerr << "\nDirect coupling calculate X vector from node invoked." << flush; exit(1); };
	virtual inline void directCoupling(bool& /*distributionFine*/, bool& /*distributionCoarse*/, Control& /*bc*/) { cerr << "\nDirect coupling from node invoked." << flush; exit(1); };
	virtual inline void interpolTimeDirectCoupling(bool& /*distributionFine*/, bool& /*distributonCoarseFuture*/, Control& /*bc*/) { cerr << "\nTime interpolation with direct coupling from node invoked." << flush; exit(1); }; //!virtual declaration for wallNode
	virtual inline void interpolTimeScalingCoarseToFine(bool& /*distributionFine*/, bool& /*distributonCoarseFuture*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nTime interpolation and scaling coarse to fine from node invoked." << flush; exit(1); }; //!virtual declaration for wallNode
	virtual inline void interpolTimeScalingCoarseToFineStrainRates(int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nTime interpolation and scaling coarse to fine from node invoked." << flush; exit(1); }; //!virtual declaration for wallNode
	virtual inline void scalingCoarseToFine(bool& /*distributionFine*/, bool& /*distributionCoarse*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nScaling coarse to fine from node invoked." << flush; exit(1); };
	virtual inline void scalingFineToCoarse(bool&/*distributionCoarse*/, bool& /*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nScaling fine to coarse from node invoked." << flush; exit(1); };
	virtual inline void interpolTimeScalingCoarseToFineMRT_GSbasis(bool& /*distributionFine*/, bool& /*distributonCoarseFuture*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nMRT time interpolation and scaling coarse to fine from node invoked." << flush; exit(1); }; //!virtual declaration for wallNode
	virtual inline void scalingCoarseToFineMRT_GSbasis(bool& /*distributionFine*/, bool& /*distributionCoarse*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nMRT scaling coarse to fine from node invoked." << flush; exit(1); };
	virtual inline void scalingFineToCoarseMRT_GSbasis(bool&/*distributionCoarse*/, bool& /*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nMRT scaling fine to coarse from node invoked." << flush; exit(1); };
	virtual inline void scalingFineToCoarseMRT_TOU_GSbasis(bool&/*distributionCoarse*/, bool& /*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nMRT scaling fine to coarse from node invoked." << flush; exit(1); };
	virtual inline void scalingFineToCoarseMRT_LAG_GSbasis(bool&/*distributionCoarse*/, bool& /*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nMRT scaling fine to coarse from node invoked." << flush; exit(1); };
	virtual inline void scalingFineToCoarseLAG(bool&/*distributionCoarse*/, bool& /*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nLAGRAVA scaling fine to coarse from node invoked." << flush; exit(1); };
	virtual inline void scalingFineToCoarseTOU(bool&/*distributionCoarse*/, bool& /*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*level*/, Node*****& /*mesh_*/, Control& /*bc*/) { cerr << "\nTOUIL scaling fine to coarse from node invoked." << flush; exit(1); };
	virtual inline void interpolSpaceMidNode(bool&/*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*d1Interpol*/, int& /*level*/, Control& /*bc_*/, Node*****& /*mesh_*/, string& /*replacementString*/) { cerr << "\nMid-node spatial interpolation from node invoked." << flush; exit(1); };
	virtual inline void interpolSpaceCenterNode(bool&/*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*d2Interface*/, int& /*level*/, Control& /*bc_*/, Node*****& /*mesh_*/, string& /*replacementString*/) { cerr << "\nCenter-node spatial interpolation from node invoked." << flush; exit(1); };
	virtual inline void interpolSpaceMidNode(bool&/*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*d1Interpol*/, int& /*level*/, Control& /*bc_*/, Node*****& /*mesh_*/) { cerr << "\nMid-node spatial interpolation from node invoked." << flush; exit(1); };
	virtual inline void interpolSpaceCenterNode(bool&/*distributionFine*/, int& /*i*/, int& /*j*/, int& /*k*/, int& /*d2Interface*/, int& /*level*/, Control& /*bc_*/, Node*****& /*mesh_*/) { cerr << "\nCenter-node spatial interpolation from node invoked." << flush; exit(1); };
	virtual inline void interpolSpaceMidNodeStrainRates(int& /*i*/, int& /*j*/, int& /*k*/, int& /*d1Interpol*/, int& /*level*/, Control& /*bc_*/, Node*****& /*mesh_*/) { cerr << "\nMid-node spatial interpolation from node invoked." << flush; exit(1); };
	virtual inline void interpolSpaceCenterNodeStrainRates(int& /*i*/, int& /*j*/, int& /*k*/, int& /*d2Interface*/, int& /*level*/, Control& /*bc_*/, Node*****& /*mesh_*/) { cerr << "\nCenter-node spatial interpolation from node invoked." << flush; exit(1); };

private:

	//!coordinates
	double x_;
	double y_;
	double z_;

protected:

	//!level
	int level_;

	//!node indices
	int indexI_;
	int indexJ_;
	int indexK_;

	//!distributions
	double distr_[2][19];

	//!distribution array for fictitous transport step
	double distrFictitious_[19];

	//!cubic Mach correction terms
	double psi_[19];

	//!cubic Mach correction terms
	double psiPreCol_[19];

	//!cubic Mach correction terms @ t + 1/2 (second order central finite difference scheme)
	double psiHalf_[19];

	//!cubic Mach correction terms @ t + 1 (second order central finite difference scheme)
	double psiFull_[19];

	//!Pi1 stress tensor
	double Pi1Tensor_[3][3];

	//!strain-rate tensor (second order central finite difference scheme)
	double strainRateTensorFD_[3][3];

	//!strain-rate tensor @ t + 1 (second order central finite difference scheme)
	double strainRateTensorFDfull_[3][3];

	//!strain-rate tensor @ t + 1/2 (second order central finite difference scheme)
	double strainRateTensorFDhalf_[3][3];

	//!pre-collision distribution
	double distrPreCol_[19];

	//!relaxation time
	double colFreq_;

	//!eddy viscosity
	double turbVisc_;

	//!density
	double densityNode_;

	//!velocity vector
	double velocityNode_[3];

	//!time step
	double timeStep_;

	//!array of pointers to adjacent nodes --> create connectivity D3Q19
	Node* neighbors_[18];

	//!partner node
	Node* partner_;

	//!first and second normal-direction neighbors of boundary nodes for LODI BCs
	Node* normalNeighbor_;
	Node* normalNeighborSnd_;
	Node* tangentNeighbor1_[2];
	Node* tangentNeighbor2_[2];
	int normalDir_;
	bool isCorner_;

};
