#pragma once

#include "Node.h"

using namespace std;

class GridCoupling : public Node
{
public:

	//!constructor
	GridCoupling() : Node() { tag_ = "gridcoupling"; ID_ = 2; }

	GridCoupling(double x, double y, double z, int indexI, int indexJ, int indexK, int level, string gridcouplingTag) : Node(x, y, z, indexI, indexJ, indexK, level)
	{
		tag_ = "gridcoupling";
		ID_ = 2;
		gridCouplingTag_ = gridcouplingTag;
	}

	//!destructor
	virtual ~GridCoupling()
	{

	}

	//!members
	virtual void interpolTimeDirectCoupling(bool& distrFi, bool& distrCoFut, Control& bc);
	virtual void directCoupling(bool& distrFi, bool& distrCo, Control& bc);
	virtual void interpolTimeScalingCoarseToFine(bool& distrFi, bool& distrCoFut, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc); //!scaling and linear time interpolation
	virtual void interpolTimeScalingCoarseToFineStrainRates(int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc); //!scaling and linear time interpolation
	virtual void interpolTimeScalingCoarseToFineMRT_GSbasis(bool& distrFi, bool& distrCoFut, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc); //!scaling and linear time interpolation
	virtual void scalingCoarseToFine(bool& distrFi, bool& distrCo, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void scalingFineToCoarse(bool& distrCo, bool& distrFi, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void scalingCoarseToFineMRT_GSbasis(bool& distrFi, bool& distrCo, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void scalingFineToCoarseMRT_GSbasis(bool& distrCo, bool& distrFi, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void scalingFineToCoarseMRT_TOU_GSbasis(bool& distrCo, bool& distrFi, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void scalingFineToCoarseMRT_LAG_GSbasis(bool& distrCo, bool& distrFi, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void scalingFineToCoarseLAG(bool& distrCo, bool& distrFi, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void scalingFineToCoarseTOU(bool& distrCo, bool& distrFi, int& i, int& j, int& k, int& level, Node*****& mesh_, Control& bc);
	virtual void interpolSpaceMidNode(bool& distrFi, int& i, int& j, int& k, int& d1Interpol, int& level, Control& bc_, Node*****& mesh_);
	virtual void interpolSpaceCenterNode(bool& distrFi, int& i, int& j, int& k, int& d2Interpol, int& level, Control& bc_, Node*****& mesh_);
	virtual void interpolSpaceMidNodeStrainRates(int& i, int& j, int& k, int& d1Interpol, int& level, Control& bc_, Node*****& mesh_);
	virtual void interpolSpaceCenterNodeStrainRates(int& i, int& j, int& k, int& d2Interpol, int& level, Control& bc_, Node*****& mesh_);
	virtual void directCouplingCalcDeltaX(double* X, double* distr, double& rho, double* velo, string& interfaceType, Control& bc);
	virtual string getGridCouplingTag(){return gridCouplingTag_;}

	virtual inline void transport(bool& distr) //!transport step (push step)
	{
		for (int dir = 0; dir < 18; dir++)
		{
            if (neighbors_[dir] != NULL)
			{
                neighbors_[dir]->setDistribution(distr, dir, distr_[!distr][dir]);
            }
		}
		distr_[distr][18] = distr_[!distr][18];
	}

	virtual inline void transportPull(bool& distr) //!transport step (pull step)
	{
		int counterDir;
		int zeroDir = 18;
		bool negDistr = !distr;

		distrFictitious_[18] = 0.;
		for (int dir = 0; dir < 18; dir++)
		{
			switch (dir)
			{
			case 0:  counterDir = 2;  break;
			case 1:  counterDir = 3;  break;
			case 2:  counterDir = 0;  break;
			case 3:  counterDir = 1;  break;
			case 4:  counterDir = 5;  break;
			case 5:  counterDir = 4;  break;
			case 6:  counterDir = 11; break;
			case 7:  counterDir = 10; break;
			case 8:  counterDir = 13; break;
			case 9:  counterDir = 12; break;
			case 10: counterDir = 7;  break;
			case 11: counterDir = 6;  break;
			case 12: counterDir = 9;  break;
			case 13: counterDir = 8;  break;
			case 14: counterDir = 15; break;
			case 15: counterDir = 14; break;
			case 16: counterDir = 17; break;
			case 17: counterDir = 16; break;
			}
			distrFictitious_[dir] = neighbors_[counterDir]->getDistribution(distr, dir);
			distrFictitious_[18] += neighbors_[dir]->getDistribution(negDistr, zeroDir);
		}
		distrFictitious_[18] /= 18.;
	}


private:

protected:

	string gridCouplingTag_;

    //!base points for compact interpolation
	Node* partnerCenterNodeAIOCTA_[8];


};
