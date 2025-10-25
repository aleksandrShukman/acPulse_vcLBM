#pragma once

#include "Node.h"

using namespace std;

class BoundaryNode : public Node
{
public:
	//!constructor
	BoundaryNode() : Node() { tag_ = "boundary"; ID_ = 1; }

	BoundaryNode(double x, double y, double z, int indexI, int indexJ, int indexK, int level) : Node(x, y, z, indexI, indexJ, indexK, level)
	{
		tag_ = "boundary";
		ID_ = 1;
	}

	//!destructor
	virtual ~BoundaryNode() {}

	virtual inline void transport(bool& alternating) //!transport step (push step)
	{
		for (int dir = 0; dir < 18; dir++)
		{
			if (neighbors_[dir] != NULL)
			{
				neighbors_[dir]->setDistribution(alternating, dir, distr_[!alternating][dir]);
			}
		}
		distr_[alternating][18] = distr_[!alternating][18];
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

};
