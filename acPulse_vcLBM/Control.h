#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

class Control
{
	public:

	//!constructor
	Control() {}

	Control(string filename);

	//!destructor
	~Control() {}

	//!members: read functions for parameters from control file
	double readValue(string filename,string keyword); //!filename: name of control file; keyword: string to be read from file; return value: (double) to corresponding keyword
	string readString(string filename,string keyword); //!filename: name of control file; keyword: string to be read from file; return value: (string) to corresponding keyword
	string readFileName(string filename,string keyword); //!filename: name of control file; keyword: string to be read from file; return value: (string) filname (a string including posibility for whitespaces)
	void readVector(string filename,string keyword,double* vector);	//!filename: name of control file; keyword: string to be read from file; return value a vector (3D!) to corresponding keyword --> the 3 components are written to array

	//!get and set members
	inline int			  getNumberOfNodesX() { return numberOfNodesX_; }
	inline int		      getNumberOfNodesY() { return numberOfNodesY_; }
	inline int			  getNumberOfNodesZ() { return numberOfNodesZ_; }
	inline int            getLevelMax() { return levelMax_; }
	inline int			  getTimeStepMax() { return timeStepMax_; }
	inline int			  getWriteInterval() { return writeInterval_; }
	inline int			  getNumberOfThreads() { return numberOfThreads_; }
	inline int		      getRefinementLayer() { return refLayer_; }
	inline int		      getRefinementLayerX() { return refLayerX_; }
	inline int            getOrderOfEquilibrium() { return equilOrder_; }
    inline void           setTime_(double time) { time_ = time; }
	inline double		  getSpacing() {return spacing_;}
	inline double		  getHybridParameter() { return hybridPar_; }
	inline double		  getDensity() {return density_;}
	inline double		  getMachNumber() {return machNumber_;}
	inline double		  getMolecularVelocity() {return molecularVelocity_;}
	inline double		  getReynoldsNumber() { return reynoldsNumber_; }
	inline double		  getSpeedOfSound() {return speedOfSound_;}
	inline double		  getKinematicViscosity() {return kinematicViscosity_;}
	inline double         getTimeStep() {return timeStep_;}
	inline double		  getTime_() { return time_; }
	inline double		  getChannelSizeX_() {return channelSizeX_;}
	inline double		  getChannelSizeY_() {return channelSizeY_;}
	inline double		  getChannelSizeZ_() {return channelSizeZ_;}
    inline double         getXsi(int& dir, int& comp) { return xsi_[dir][comp]; } //!unsafe: reference to local variable is passed; watch out when using
    inline double          getXsi(unsigned short& dir, unsigned short& comp) { return xsi_[dir][comp]; } //!unsafe: reference to local variable is passed; watch out when using
	inline double*        getXsiX() {return xsiX_;} //!unsafe: reference to local variable is passed; watch out when using
	inline double*        getXsiY() {return xsiY_;} //!unsafe: reference to local variable is passed; watch out when using
	inline double*        getXsiZ() {return xsiZ_;} //!unsafe: reference to local variable is passed; watch out when using
	inline double*		  getXsi(int& dir) { return xsi_[dir]; } //!unsafe: reference to local variable is passed; watch out when using
	inline double*		  getXsi(unsigned short& dir) { return xsi_[dir]; } //!unsafe: reference to local variable is passed; watch out when using
	inline double*        getWeightFactors() {return weightFactors_;} //!unsafe: reference to loacl variable is passed; watch out when using
	inline double		  getLODIK1relaxation() { return lodiK1relaxation_; }
	inline double		  getLODIK2relaxation() { return lodiK2relaxation_; }
    virtual inline double getCollisionFrequency() { return collisionFrequency_; }
    virtual inline double getOverlap() { return overlap_; }
	inline string         getCaseName() { return caseName_; }
	inline string		  getCollisionModel() { return collisionModel_; }
	inline string		  getCouplingReplacementC2F() { return couplingReplacementC2F_; }
	inline string		  getCouplingReplacementF2C() { return couplingReplacementF2C_; }
    inline string		  getHangingNodeReconstruction() { return hangingNodeReconstruction_; }
	inline string		  getScalingStyle() { return scalingStyle_; }
	inline string		  getReflectionHandling() { return reflecting_; }
	inline string		  getTransverseInclusion() { return transverse_; }
	inline string		  getCubicMachCorrection() { return cubicMachCorrection_; }
	inline string		  getF2Cfilter() { return F2Cfilter_; }
    inline string		  getCouplingProcedure() { return couplingProcedure_; }

	private:

	//!grid variables
	int    levelMax_; //!max. level of mesh
	int    numberOfNodesX_, numberOfNodesY_, numberOfNodesZ_; //!number of nodes in x, y and z direction
	int    numberOfCellsX_, numberOfCellsY_, numberOfCellsZ_; //!number of cells in x, y and z direction
	int	   refLayer_; //!number of refined coarse cells at YZ boundaries
	int	   refLayerX_; //!number of refined coarse cells at inlet
	int    totalNumberOfNodes_;
	int	   equilOrder_;
	double spacing_; //!mesh spacing in level 0
	double channelSizeX_, channelSizeY_, channelSizeZ_; //!channel dimensions in x, y and z direction
	double overlap_; //!number of coarse overlap-cells

	//!hybrid recursive-regularization
	double hybridPar_; //!hybridization parameter

	//!lattice variables
	double xsiX_[19], xsiY_[19], xsiZ_[19]; //!components of molecular velocity vetor (D3Q19 lattice)
	double xsi_[19][3]; //molecular velocity vetor (D3Q19 lattice)
	double weightFactors_[19]; //!lattice weight factors (D3Q19)

	//!physical variables
	double density_;				//!fluid density
	double externalAcceleration_[3];//!external acceleration
	double machNumber_;				//!global Mach number
	double molecularVelocity_;		//!molecular velocity (delta_x/delta_t)
	double speedOfSound_;			//!speed of sound
	double kinematicViscosity_;		//!kinematic velocity
	double timeStep_;				//!time step in level 0
	double collisionFrequency_;	    //!dimensionless collision frequency on level 0
	double reynoldsNumber_;		    //!reynolds number
	double lodiK1relaxation_;		//!LODI K1 relaxation factor
	double lodiK2relaxation_;		//!LODI K2 relaxation factor

	//!control parameters
	int    timeStepMax_; //!number of time steps to run(in level 0)
	int    writeInterval_; //!number of time steps between to consecutive result files
	int    numberOfThreads_; //!number of cpus for SMP
    double time_; //!physical time
	string caseName_; //!absolute path to the file to be saved
	string collisionModel_;
	string scalingStyle_;
    string couplingReplacementC2F_;
	string couplingReplacementF2C_;
    string hangingNodeReconstruction_;
	string forceModel_;
	string reflecting_;
	string transverse_;
	string cubicMachCorrection_;
	string F2Cfilter_;
    string couplingProcedure_;
};
