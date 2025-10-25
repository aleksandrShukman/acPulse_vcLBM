#pragma once

#include "Node.h"
#include "FluidNode.h"
#include "BoundaryNode.h"
#include "WallNode.h"
#include "Control.h"
#include "GridCoupling.h"
#include <omp.h>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

class Lattice
{

public:

	//!constructors

	Lattice(Control& bc) : bc_(bc) //!copy reference
	{
		//!create mesh data structure
		int maxLevel = bc.getLevelMax(); //!maximum level
		int numberOfLevels = maxLevel + 1; //!number of levels

		mesh_ = new Node****[numberOfLevels]; //!create structure for every refinement level

		for (int level = 0; level <= maxLevel; level++)
		{
			int imax = (bc.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
			int jmax = (bc.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
			int kmax = (bc.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

			int numberOfNodesX = imax + 1; //!total number of nodes in x direction in ->this level
			int numberOfNodesY = jmax + 1; //!total number of nodes in y direction in ->this level
			int numberOfNodesZ = kmax + 1; //!total number of nodes in z direction in ->this level

			mesh_[level] = new Node***[numberOfNodesX]; //!create arrays in x direction for every level

			for(int i = 0; i <= imax; i++)
			{
				mesh_[level][i] = new Node**[numberOfNodesY]; //!create arrays in y direction for every level

				for(int j = 0; j <= jmax; j++)
				{
					mesh_[level][i][j] = new Node*[numberOfNodesZ]; //!create arrays in z direction for every level
				} //!end j
			} //!end i
		} //!end level

		createLattice(); //!create nodes
		setConnectivity(); //!define connectivity

	} //!end constructor


	//!destructor
	~Lattice() //!release all dynamically allocated memory
	{
		if (mesh_ != NULL)
		{
			int maxLevel = bc_.getLevelMax(); //!maximum level
			int numberOfLevels = maxLevel + 1; //!number of levels

			for (int level = 0; level <= maxLevel; level++)
			{
				if (mesh_[level])
				{
					int imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
					int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
					int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

					for(int i = 0; i <= imax; i++)
					{
						if (mesh_[level][i])
						{
							for(int j = 0; j <= jmax; j++)
							{
								if (mesh_[level][i][j]) 
								{
									for(int k = 0; k <= kmax; k++)
									{
										if (mesh_[level][i][j][k]) delete mesh_[level][i][j][k];
									}
									delete[] mesh_[level][i][j]; //!delete arrays in z direction for every level
								}
							} //!end j
							delete[] mesh_[level][i]; //!delete arrays in y direction for every level
						}
					} //!end i
					delete[] mesh_[level]; //!delete arrays in x direction for every level
				}
			} //!end level
			delete[] mesh_; //!delete levels
		} 
	} //!end denstructor

	//!member
	virtual void createLattice();
	virtual void setConnectivity();
	virtual void solve(bool &alternating);
	virtual void solve(bool& alternatingL0, bool& alternatingL1);
	virtual void timeSteppingSRT(int& level, bool &alternating);
	virtual void timeSteppingMRT_GSbasis(int& level, bool& alternating);
	virtual void timeSteppingMRT_RMbasis(int& level, bool& alternating);
	virtual void timeSteppingRR(int& level, bool& alternating);
	virtual void timeSteppingHRR(int& level, bool& alternating);
	virtual void nestedTimeSteppingSRT(bool &alternatingL1,bool &alternatingL0, Control &bc_);
	virtual void nestedTimeSteppingMRT_GSbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void nestedTimeSteppingMRT_RMbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void nestedTimeSteppingRR(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void nestedTimeSteppingHRR(bool& alternatingL1, bool& alternatingL0, Control& bc_);
	virtual void writeResultsVTK(string filename, bool &alternating, bool& alternatingL1);
	virtual void initialization(Control& bc);
	
	virtual inline double getDistribution(bool& alternating, int& dir, int& level, int& i, int& j, int& k) { if (dir < 19) return mesh_[level][i][j][k]->getDistribution(alternating, dir); else return NULL; };
	virtual double readValue(string filename, string keyword); //!filename: name of initialization file; keyword: string to be read from file; return value: (double) to corresponding keyword

private:
	Control &bc_; //!reference to instance of control class
	Node***** mesh_; //!pointer for creating mesh data structure
};
