#pragma once

#include "Node.h"

using namespace std;

class WallNode: public Node
{

public:

	//!constructor
	WallNode() : Node() { tag_ = "wall"; ID_ = 3; } //!standard

	WallNode(double x, double y, double z, long long indexI, long long indexJ, long long indexK, long long level) : Node(x, y, z, indexI, indexJ, indexK, level)
	{
		tag_ = "wall";
		ID_ = 3;
	}

	//!destructor
	virtual ~WallNode() 
	{

	}

	virtual inline void halfWayBounceBack(bool& alternating, Control& bc)
	{
		for (int dir = 0; dir < 18; dir++) 
		{
			if (neighbors_[dir] != NULL)
			{
				//!calc opposite direction for bounce-back
				int counterDir = 18;
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

				//!perform bounce-back
				neighbors_[dir]->setDistribution(alternating, dir, distr_[alternating][counterDir]);
			}
		}
	}//!end bounceBack


	virtual inline void frictionless1stOrder(bool& alternating, Control& bc) //works only if boundary is plane and xz plane (i.e. plane has a y-normal)
	{
		//!calc neighbor direction in y-direction of fluid neighbor
		int wallDirection = 1; //!fluid neighbor is north of wall node
		if (this->getIndexJ() > 0) wallDirection = 3; //!fluid neighbor is south of wall node

		for (int dir = 0; dir < 18; dir++)
		{
			if (neighbors_[dir] != NULL)
			{
				//!calc outgoing direction for specular reflection bounce back
				//!indices: 0:E, 1:N, 2:W, 3:S, 4:T, 5:B, 6:TE, 7:BE, 8:TN, 9:BN, 10: TW, 11:BW, 12:TS, 13:BS, 14:NW, 15:SE, 16:NE, 17: SW, 18: Rest (==D3Q19)
				//!consider only walls with y-normal
				int counterDir = 19;
				int counterDirFrictionless = 19;
				//!determine counterDir (i.e. incoming direction of distribution)
				switch (dir)
				{
				case 0: counterDir = 2; break; //!east to west
				case 1: counterDir = 3; break; //!north to south
				case 2: counterDir = 0; break; //!west to east
				case 3: counterDir = 1; break; //!south to north
				case 4: counterDir = 5; break; //!top to bottom
				case 5: counterDir = 4; break; //!bottom to top
				case 6: counterDir = 11; break; //!top east to bottom west
				case 7: counterDir = 10; break; //!bottom east to top west
				case 8: counterDir = 13; break; //!top north to bottom south
				case 9: counterDir = 12; break; //!bottom north to top south
				case 10: counterDir = 7; break; //!top west to bottom east
				case 11: counterDir = 6; break; //!bottom west to top east
				case 12: counterDir = 9; break; //!top south to bottom north
				case 13: counterDir = 8; break; //!bottom south to top north
				case 14: counterDir = 15; break; //!north west to south east
				case 15: counterDir = 14; break; //!south east to north west
				case 16: counterDir = 17; break; //!north east to south west
				case 17: counterDir = 16; break; //!south west to north east
				}
				//!determine outgoing direction
				switch (counterDir)
				{
				case 1: counterDirFrictionless = 3; break; //!north to south
				case 3: counterDirFrictionless = 1; break; //!south to north
				case 8: counterDirFrictionless = 12; break; //!top north to top south
				case 9: counterDirFrictionless = 13; break; //!bottom north to bottom south
				case 12: counterDirFrictionless = 8; break; //!top south to top north
				case 13: counterDirFrictionless = 9; break; //!bottom south to bottom north
				case 14: counterDirFrictionless = 17; break; //!north west to south west
				case 15: counterDirFrictionless = 16; break; //!south east to north east
				case 16: counterDirFrictionless = 15; break; //!north east to south east
				case 17: counterDirFrictionless = 14; break; //!south west to north west
				}

				//!reflect normal component, remain tangential component
				double temp = this->getDistribution(alternating, counterDir);

				neighbors_[wallDirection]->setDistribution(alternating, counterDirFrictionless, temp);
				distr_[alternating][counterDirFrictionless] = temp;
			}
		}
	} //!end frictionless1stOrder

private:

};
