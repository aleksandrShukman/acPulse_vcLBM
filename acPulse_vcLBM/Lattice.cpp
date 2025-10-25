#define _USE_MATH_DEFINES

#include "Lattice.h"
#include "math.h"

void Lattice::createLattice()
{
	//!create nodes in mesh structure
	//!mesh topology: bottom and top layer (+/- z normal) of mesh are walls (channel walls)
	//!left and right boundaries (+/- y normal) are periodic boundaries
	//!inlet and outlet boundaries (+/- x normal) are periodic boundaries
	//!remaining nodes are fluid

	int refLayer, refLayerX, overlap;

	//!create nodes in array
	for (int level = 0; level <= bc_.getLevelMax(); level++)
	{
		int imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
		int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
		int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level

		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		overlap = bc_.getOverlap();

		for (int k = 0; k <= kmax; k++)
		{
			for (int j = 0; j <= jmax; j++)
			{
				for (int i = 0; i <= imax; i++)
				{
					double xCoord = i * (bc_.getSpacing() / pow(2., level));
					double yCoord = j * (bc_.getSpacing() / pow(2., level));
					double zCoord = k * (bc_.getSpacing() / pow(2., level));

					//!level 0 with refinement
					if ((level == 0) && (bc_.getLevelMax() > 0))
					{
						/*FACES*/
						//!gridcoupling-nodes in-plane bottom interface
						if ((k == refLayer - overlap) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap)) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "bottom"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane top interface
						else if ((k == kmax - refLayer + overlap) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap)) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "top"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane left interface
						else if ((j == refLayer - overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap)) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane  right interface
						else if ((j == jmax - refLayer + overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap)) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane front interface
						else if ((i == refLayerX - overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap)) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane back interface
						else if ((i == imax - refLayer + overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap)) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back"); //!if node is at interface --> gridcoupling node


						/*EDGES*/
						//!gridcoupling-nodes bottom left edge interface
						else if ((k == refLayer - overlap) && (j == refLayer - overlap) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "bottom_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes bottom right edge interface
						else if ((k == refLayer - overlap) && (j == jmax - refLayer + overlap) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "bottom_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes top left edge interface
						else if ((k == kmax - refLayer + overlap) && (j == refLayer - overlap) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "top_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes top right edge interface
						else if ((k == kmax - refLayer + overlap) && (j == jmax - refLayer + overlap) && ((i > refLayerX - overlap) && (i < imax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "top_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front bottom edge interface
						else if ((i == refLayerX - overlap) && (k == refLayer - overlap) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_bottom"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front top edge interface
						else if ((i == refLayerX - overlap) && (k == kmax - refLayer + overlap) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_top"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back bottom edge interface
						else if ((i == imax - refLayer + overlap) && (k == refLayer - overlap) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_bottom"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back top edge interface
						else if ((i == imax - refLayer + overlap) && (k == kmax - refLayer + overlap) && ((j > refLayer - overlap) && (j < jmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_top"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front left edge interface
						else if ((i == refLayerX - overlap) && (j == refLayer - overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front right edge interface
						else if ((i == refLayerX - overlap) && (j == jmax - refLayer + overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back left edge interface
						else if ((i == imax - refLayer + overlap) && (j == refLayer - overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back right edge interface
						else if ((i == imax - refLayer + overlap) && (j == jmax - refLayer + overlap) && ((k > refLayer - overlap) && (k < kmax - refLayer + overlap))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_right"); //!if node is at interface --> gridcoupling node


						/*CORNERS*/
						//!gridcoupling-nodes front bottom left corner interface
						else if ((i == refLayerX - overlap) && (k == refLayer - overlap) && (j == refLayer - overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_bottom_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front bottom right corner interface
						else if ((i == refLayerX - overlap) && (k == refLayer - overlap) && (j == jmax - refLayer + overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_bottom_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front top left corner interface
						else if ((i == refLayerX - overlap) && (k == kmax - refLayer + overlap) && (j == refLayer - overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_top_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front top right corner interface
						else if ((i == refLayerX - overlap) && (k == kmax - refLayer + overlap) && (j == jmax - refLayer + overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_top_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back bottom left corner interface
						else if ((i == imax - refLayer + overlap) && (k == refLayer - overlap) && (j == refLayer - overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_bottom_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back bottom right corner interface
						else if ((i == imax - refLayer + overlap) && (k == refLayer - overlap) && (j == jmax - refLayer + overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_bottom_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back top left corner interface
						else if ((i == imax - refLayer + overlap) && (k == kmax - refLayer + overlap) && (j == refLayer - overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_top_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back top right corner interface
						else if ((i == imax - refLayer + overlap) && (k == kmax - refLayer + overlap) && (j == jmax - refLayer + overlap)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_top_right"); //!if node is at interface --> gridcoupling node


						//!standard fluid nodes
						else mesh_[level][i][j][k] = new FluidNode(xCoord, yCoord, zCoord, i, j, k, level); //!all other nodes are fluid
					}

					//!level 0 without refinement
					else if (bc_.getLevelMax() == 0)
					{
						//!boundary nodes
						if (i == 0 || i == imax || j == 0 || j == jmax || k == 0 || k == kmax) mesh_[level][i][j][k] = new BoundaryNode(xCoord, yCoord, zCoord, i, j, k, level); //!if node is at front, back, bottom, top, left or right side --> boundary node

						//!standard fluid nodes
						else mesh_[level][i][j][k] = new FluidNode(xCoord, yCoord, zCoord, i, j, k, level); //!all other nodes are fluid
					}

					//!level 1
					else if (level == 1)
					{
						//!boundary nodes
						if (i == 0 || i == imax || j == 0 || j == jmax || k == 0 || k == kmax) mesh_[level][i][j][k] = new BoundaryNode(xCoord, yCoord, zCoord, i, j, k, level); //!if node is at front, back, bottom, top, left or right side --> boundary node

						/*FACES*/
						//!gridcoupling-nodes in-plane bottom interface
						else if ((k == refLayer) && ((j > refLayer) && (j < jmax - refLayer)) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "bottom"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane top interface
						else if ((k == kmax - refLayer) && ((j > refLayer) && (j < jmax - refLayer)) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "top"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane left interface
						else if ((j == refLayer) && ((k > refLayer) && (k < kmax - refLayer)) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane  right interface
						else if ((j == jmax - refLayer) && ((k > refLayer) && (k < kmax - refLayer)) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane front interface
						else if ((i == refLayerX) && ((k > refLayer) && (k < kmax - refLayer)) && ((j > refLayer) && (j < jmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes in-plane back interface
						else if ((i == imax - refLayer) && ((k > refLayer) && (k < kmax - refLayer)) && ((j > refLayer) && (j < jmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back"); //!if node is at interface --> gridcoupling node


						/*EDGES*/
						//!gridcoupling-nodes bottom left edge interface
						else if ((k == refLayer) && (j == refLayer) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "bottom_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes bottom right edge interface
						else if ((k == refLayer) && (j == jmax - refLayer) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "bottom_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes top left edge interface
						else if ((k == kmax - refLayer) && (j == refLayer) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "top_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes top right edge interface
						else if ((k == kmax - refLayer) && (j == jmax - refLayer) && ((i > refLayerX) && (i < imax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "top_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front bottom edge interface
						else if ((i == refLayerX) && (k == refLayer) && ((j > refLayer) && (j < jmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_bottom"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front top edge interface
						else if ((i == refLayerX) && (k == kmax - refLayer) && ((j > refLayer) && (j < jmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_top"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back bottom edge interface
						else if ((i == imax - refLayer) && (k == refLayer) && ((j > refLayer) && (j < jmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_bottom"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back top edge interface
						else if ((i == imax - refLayer) && (k == kmax - refLayer) && ((j > refLayer) && (j < jmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_top"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front left edge interface
						else if ((i == refLayerX) && (j == refLayer) && ((k > refLayer) && (k < kmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front right edge interface
						else if ((i == refLayerX) && (j == jmax - refLayer) && ((k > refLayer) && (k < kmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back left edge interface
						else if ((i == imax - refLayer) && (j == refLayer) && ((k > refLayer) && (k < kmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back right edge interface
						else if ((i == imax - refLayer) && (j == jmax - refLayer) && ((k > refLayer) && (k < kmax - refLayer))) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_right"); //!if node is at interface --> gridcoupling node


						/*CORNERS*/
						//!gridcoupling-nodes front bottom left corner interface
						else if ((i == refLayerX) && (k == refLayer) && (j == refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_bottom_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front bottom right corner interface
						else if ((i == refLayerX) && (k == refLayer) && (j == jmax - refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_bottom_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front top left corner interface
						else if ((i == refLayerX) && (k == kmax - refLayer) && (j == refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_top_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes front top right corner interface
						else if ((i == refLayerX) && (k == kmax - refLayer) && (j == jmax - refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "front_top_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back bottom left corner interface
						else if ((i == imax - refLayer) && (k == refLayer) && (j == refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_bottom_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back bottom right corner interface
						else if ((i == imax - refLayer) && (k == refLayer) && (j == jmax - refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_bottom_right"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back top left corner interface
						else if ((i == imax - refLayer) && (k == kmax - refLayer) && (j == refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_top_left"); //!if node is at interface --> gridcoupling node

						//!gridcoupling-nodes back top right corner interface
						else if ((i == imax - refLayer) && (k == kmax - refLayer) && (j == jmax - refLayer)) mesh_[level][i][j][k] = new GridCoupling(xCoord, yCoord, zCoord, i, j, k, level, "back_top_right"); //!if node is at interface --> gridcoupling node

						//!standard fluid nodes
						else mesh_[level][i][j][k] = new FluidNode(xCoord, yCoord, zCoord, i, j, k, level);	//!all other nodes are fluid
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end level
} //!end createLattice


void Lattice::setConnectivity()
{
	//!create connectivity, i.e., set neighbor pointer of node to the according neighbor node
	//!for wall nodes the neighbor pointer is set to NULL if neighbor is missing
	//!neighbor pointer for periodic nodes are set to periodic neighbor to achieve periodicity condition

	//!set connectivity for fluid nodes
	for (int level = 0; level <= bc_.getLevelMax(); level++)
	{
		int imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
		int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
		int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level

		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction except inlet and outlet boundaries
		{
			for (int j = 1; j <= jmax - 1; j++) //!all nodes except left and right wall
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except bottom and top wall
				{
					for (int dir = 0; dir < 18; dir++)
					{
						int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
                        switch (dir) //!calc indices of neighbor
						{
						case 0:
							if (mesh_[level][i][j][k]->getTag() == "gridcoupling")
							{
								if (((mesh_[level][i][j][k]->getGridCouplingTag() == "front") && (level == 1)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top_right") && (level == 0)))
								{
									Node* neighbor = NULL;
									mesh_[level][i][j][k]->setNeighbor(dir, neighbor);
									continue;
								}
								else
								{
									iNeighbor = i + 1;
									jNeighbor = j;
									kNeighbor = k;
									break;
								}
							}
							else
							{
								iNeighbor = i + 1;
								jNeighbor = j;
								kNeighbor = k;
								break;
							}
						case 1:
							if (mesh_[level][i][j][k]->getTag() == "gridcoupling")
							{
								if (((mesh_[level][i][j][k]->getGridCouplingTag() == "left") && (level == 1)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "bottom_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "top_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top_right") && (level == 0)))
								{
									Node* neighbor = NULL;
									mesh_[level][i][j][k]->setNeighbor(dir, neighbor);
									continue;
								}
								else
								{
									iNeighbor = i;
									jNeighbor = j + 1;
									kNeighbor = k;
									break;
								}
							}
							else
							{
								iNeighbor = i;
								jNeighbor = j + 1;
								kNeighbor = k;
								break;
							}
						case 2:
							if (mesh_[level][i][j][k]->getTag() == "gridcoupling")
							{
								if (((mesh_[level][i][j][k]->getGridCouplingTag() == "back") && (level == 1)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top_right") && (level == 0)))
								{
									Node* neighbor = NULL;
									mesh_[level][i][j][k]->setNeighbor(dir, neighbor);
									continue;
								}
								else
								{
									iNeighbor = i - 1;
									jNeighbor = j;
									kNeighbor = k;
									break;
								}
							}
							else
							{
								iNeighbor = i - 1;
								jNeighbor = j;
								kNeighbor = k;
								break;
							}
						case 3:
							if (mesh_[level][i][j][k]->getTag() == "gridcoupling")
							{
								if (((mesh_[level][i][j][k]->getGridCouplingTag() == "right") && (level == 1)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "top_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top_left") && (level == 0)))
								{
									Node* neighbor = NULL;
									mesh_[level][i][j][k]->setNeighbor(dir, neighbor);
									continue;
								}
								else
								{
									iNeighbor = i;
									jNeighbor = j - 1;
									kNeighbor = k;
									break;
								}
							}
							else
							{
								iNeighbor = i;
								jNeighbor = j - 1;
								kNeighbor = k;
								break;
							}
						case 4:
							if (mesh_[level][i][j][k]->getTag() == "gridcoupling")
							{
								if (((mesh_[level][i][j][k]->getGridCouplingTag() == "bottom") && (level == 1)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "top") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "top_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "top_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_top_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_top_right") && (level == 0)))
								{
									Node* neighbor = NULL;
									mesh_[level][i][j][k]->setNeighbor(dir, neighbor);
									continue;
								}
								else
								{
									iNeighbor = i;
									jNeighbor = j;
									kNeighbor = k + 1;
									break;
								}
							}
							else
							{
								iNeighbor = i;
								jNeighbor = j;
								kNeighbor = k + 1;
								break;
							}
						case 5:
							if (mesh_[level][i][j][k]->getTag() == "gridcoupling")
							{
								if (((mesh_[level][i][j][k]->getGridCouplingTag() == "top") && (level == 1)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "bottom") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "bottom_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "front_bottom_right") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom_left") && (level == 0)) || ((mesh_[level][i][j][k]->getGridCouplingTag() == "back_bottom_right") && (level == 0)))
								{
									Node* neighbor = NULL;
									mesh_[level][i][j][k]->setNeighbor(dir, neighbor);
									continue;
								}
								else
								{
									iNeighbor = i;
									jNeighbor = j;
									kNeighbor = k - 1;
									break;
								}
							}
							else
							{
								iNeighbor = i;
								jNeighbor = j;
								kNeighbor = k - 1;
								break;
							}
						case 6:
							iNeighbor = i + 1;
							jNeighbor = j;
							kNeighbor = k + 1;
							break;
						case 7:
							iNeighbor = i + 1;
							jNeighbor = j;
							kNeighbor = k - 1;
							break;
						case 8:
							iNeighbor = i;
							jNeighbor = j + 1;
							kNeighbor = k + 1;
							break;
						case 9:
							iNeighbor = i;
							jNeighbor = j + 1;
							kNeighbor = k - 1;
							break;
						case 10:
							iNeighbor = i - 1;
							jNeighbor = j;
							kNeighbor = k + 1;
							break;
						case 11:
							iNeighbor = i - 1;
							jNeighbor = j;
							kNeighbor = k - 1;
							break;
						case 12:
							iNeighbor = i;
							jNeighbor = j - 1;
							kNeighbor = k + 1;
							break;
						case 13:
							iNeighbor = i;
							jNeighbor = j - 1;
							kNeighbor = k - 1;
							break;
						case 14:
							iNeighbor = i - 1;
							jNeighbor = j + 1;
							kNeighbor = k;
							break;
						case 15:
							iNeighbor = i + 1;
							jNeighbor = j - 1;
							kNeighbor = k;
							break;
						case 16:
							iNeighbor = i + 1;
							jNeighbor = j + 1;
							kNeighbor = k;
							break;
						case 17:
							iNeighbor = i - 1;
							jNeighbor = j - 1;
							kNeighbor = k;
							break;
						} //!end switch

						Node* neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						//!check if neighbor was found
						if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
						mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
					} //!end dir
					if ((mesh_[level][i][j][k]->getTag() == "gridcoupling") && (level == 0))
					{
						int n = 2;
						int ii = i * n;
						int jj = j * n;
						int kk = k * n;
						Node* partner = mesh_[level + 1][ii][jj][kk];
						mesh_[level][i][j][k]->setPartner(partner);
					}
					else if ((mesh_[level][i][j][k]->getTag() == "gridcoupling") && (level > 0))
					{
						int n = 2;
						int ii = i / n;
						int jj = j / n;
						int kk = k / n;
						Node* partner = mesh_[level - 1][ii][jj][kk];
						mesh_[level][i][j][k]->setPartner(partner);
					}
					else
					{
						Node* partner = NULL;
						mesh_[level][i][j][k]->setPartner(partner);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end level

	 //!set connectivity for boundary nodes
	for (int level = 0; level <= bc_.getLevelMax(); level++)
	{
		int imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
		int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
		int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level

		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction except inlet and outlet edges
		{
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in y direction except corners
			{
				//!do bottom boundary
				int k = 0;

				for (int dir = 0; dir < 18; dir++)
				{
					Node* neighbor = NULL;
					Node* temp = NULL;
					Node* temp2 = NULL;
					int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
					int iTNeighbor, jTNeighbor, kTNeighbor; //!node indices of next neighbor in this direction
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @ bottom only nodes in top directions matter
					{
						//!in-plane neighbors for transverse LODI
					case 0: iTNeighbor = i + 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 2: iTNeighbor = i - 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;
					case 1: iTNeighbor = i; jTNeighbor = j + 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 3: iTNeighbor = i; jTNeighbor = j - 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;

					case 4:  iNeighbor = i;     jNeighbor = j;     kNeighbor = k + 1;
						temp = mesh_[level][i][j][k + 1];
						temp2 = mesh_[level][i][j][k + 2];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						break;
					case 6:  iNeighbor = i + 1; jNeighbor = j;     kNeighbor = k + 1; break;
					case 8:	 iNeighbor = i;     jNeighbor = j + 1; kNeighbor = k + 1; break;
					case 10: iNeighbor = i - 1; jNeighbor = j;     kNeighbor = k + 1; break;
					case 12: iNeighbor = i;     jNeighbor = j - 1; kNeighbor = k + 1; break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch

					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];

					//!check if neighbor was found
					if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
					mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
				} //!end dir

				//!do top boundary
				k = kmax;

				for (int dir = 0; dir < 18; dir++)
				{
					Node* neighbor = NULL;
					Node* temp = NULL;
					Node* temp2 = NULL;
					int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
					int iTNeighbor, jTNeighbor, kTNeighbor; //!node indices of next neighbor in this direction
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @top only nodes in bottom directions matter
					{
						//!in-plane neighbors for transverse LODI
					case 0: iTNeighbor = i + 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 2: iTNeighbor = i - 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;
					case 1: iTNeighbor = i; jTNeighbor = j + 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 3: iTNeighbor = i; jTNeighbor = j - 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;

					case 5:  iNeighbor = i;     jNeighbor = j;     kNeighbor = k - 1;
						temp = mesh_[level][i][j][k - 1];
						temp2 = mesh_[level][i][j][k - 2];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						break;
					case 7:  iNeighbor = i + 1; jNeighbor = j;     kNeighbor = k - 1; break;
					case 9:  iNeighbor = i;     jNeighbor = j + 1; kNeighbor = k - 1; break;
					case 11: iNeighbor = i - 1; jNeighbor = j;     kNeighbor = k - 1; break;
					case 13: iNeighbor = i;     jNeighbor = j - 1; kNeighbor = k - 1; break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch

					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];

					//!check if neighbor was found
					if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
					mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
				} //!end dir
			} //!end j

			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except corners
			{
				//!do left boundary
				int j = 0;

				for (int dir = 0; dir < 18; dir++)
				{
					Node* neighbor = NULL;
					Node* temp = NULL;
					Node* temp2 = NULL;
					int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
					int iTNeighbor, jTNeighbor, kTNeighbor; //!node indices of next neighbor in this direction
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @left only nodes in right directions matter
					{
						//!in-plane neighbors for transverse LODI
					case 0: iTNeighbor = i + 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 2: iTNeighbor = i - 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;
					case 4: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k + 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 5: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k - 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;

					case 1:  iNeighbor = i;		jNeighbor = j + 1;	kNeighbor = k;
						temp = mesh_[level][i][j + 1][k];
						temp2 = mesh_[level][i][j + 2][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						break;
					case 16: iNeighbor = i + 1; jNeighbor = j + 1;	kNeighbor = k;     break;
					case 8:	 iNeighbor = i;		jNeighbor = j + 1;	kNeighbor = k + 1; break;
					case 14: iNeighbor = i - 1; jNeighbor = j + 1;	kNeighbor = k;     break;
					case 9:  iNeighbor = i;     jNeighbor = j + 1;  kNeighbor = k - 1; break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch

					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];

					//!check if neighbor was found
					if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
					mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
				} //!end dir

				//!do right boundary
				j = jmax;

				for (int dir = 0; dir < 18; dir++)
				{
					Node* neighbor = NULL;
					Node* temp = NULL;
					Node* temp2 = NULL;
					int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
					int iTNeighbor, jTNeighbor, kTNeighbor; //!node indices of next neighbor in this direction
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @right only nodes in left directions matter
					{
						//!in-plane neighbors for transverse LODI
					case 0: iTNeighbor = i + 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 2: iTNeighbor = i - 1; jTNeighbor = j; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;
					case 4: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k + 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 5: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k - 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;

					case 3:  iNeighbor = i;     jNeighbor = j - 1;	kNeighbor = k;
						temp = mesh_[level][i][j - 1][k];
						temp2 = mesh_[level][i][j - 2][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						break;
					case 15: iNeighbor = i + 1; jNeighbor = j - 1;	kNeighbor = k;     break;
					case 12: iNeighbor = i;		jNeighbor = j - 1;	kNeighbor = k + 1; break;
					case 17: iNeighbor = i - 1; jNeighbor = j - 1;  kNeighbor = k;     break;
					case 13: iNeighbor = i;		jNeighbor = j - 1;	kNeighbor = k - 1; break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch

					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];

					//!check if neighbor was found
					if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
					mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
				} //!end dir
			} //!end k

			//!do edge nodes
			//!bottom left edge
			int k = 0;
			int j = 0;

			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* neighborSnd = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @ bottom only nodes in top directions matter
				{
				case 1:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 8: iNeighbor = i; jNeighbor = j + 1; kNeighbor = k + 1;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor][jNeighbor + 1][kNeighbor + 1];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir

			//!top left edge
			k = kmax;
			j = 0;

			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
				{
				case 1:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 9: iNeighbor = i; jNeighbor = j + 1; kNeighbor = k - 1;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor][jNeighbor + 1][kNeighbor - 1];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir

			//!bottom right edge
			k = 0;
			j = jmax;

			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @ bottom only nodes in top directions matter
				{
				case 3:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 12: iNeighbor = i; jNeighbor = j - 1; kNeighbor = k + 1;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor][jNeighbor - 1][kNeighbor + 1];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir

			//!top right edge
			k = kmax;
			j = jmax;

			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
				{
				case 3:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 13: iNeighbor = i; jNeighbor = j - 1; kNeighbor = k - 1;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor][jNeighbor - 1][kNeighbor - 1];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end i

		//!do inlet nodes
		int i = 0;

		for (int j = 1; j <= jmax - 1; j++)
		{
			for (int k = 1; k <= kmax - 1; k++)
			{
				for (int dir = 0; dir < 18; dir++)
				{
					Node* neighbor = NULL;
					Node* temp = NULL;
					Node* temp2 = NULL;
					int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
					int iTNeighbor, jTNeighbor, kTNeighbor; //!node indices of next neighbor in this direction
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
						//!in-plane neighbors for transverse LODI
					case 1: iTNeighbor = i; jTNeighbor = j + 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 3: iTNeighbor = i; jTNeighbor = j - 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;
					case 4: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k + 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0); continue;
					case 5: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k - 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1); continue;

					case 0:  iNeighbor = i + 1; jNeighbor = j; kNeighbor = k;
						temp = mesh_[level][i + 1][j][k];
						temp2 = mesh_[level][i + 2][j][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						break;
					case 6:  iNeighbor = i + 1; jNeighbor = j;		kNeighbor = k + 1; break;
					case 7:  iNeighbor = i + 1; jNeighbor = j;		kNeighbor = k - 1; break;
					case 15: iNeighbor = i + 1; jNeighbor = j - 1;	kNeighbor = k;	   break;
					case 16: iNeighbor = i + 1; jNeighbor = j + 1;	kNeighbor = k;     break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch

					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];

					//!check if neighbor was found
					if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
					mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
				} //!end dir
			} //!end k
		} //!end j

		//!do outlet nodes
		i = imax;

		for (int j = 1; j <= jmax - 1; j++)
		{
			for (int k = 1; k <= kmax - 1; k++)
			{
				for (int dir = 0; dir < 18; dir++)
				{
					Node* neighbor = NULL;
					Node* temp = NULL;
					Node* temp2 = NULL;
					int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
					int iTNeighbor, jTNeighbor, kTNeighbor; //!node indices of next neighbor in this direction
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
						//!in-plane neighbors for transverse LODI
					case 1: iTNeighbor = i; jTNeighbor = j + 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0);  continue;
					case 3: iTNeighbor = i; jTNeighbor = j - 1; kTNeighbor = k;
						mesh_[level][i][j][k]->setTangentNeighbor1(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1);  continue;
					case 4: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k + 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 0);  continue;
					case 5: iTNeighbor = i; jTNeighbor = j; kTNeighbor = k - 1;
						mesh_[level][i][j][k]->setTangentNeighbor2(mesh_[level][iTNeighbor][jTNeighbor][kTNeighbor], 1);  continue;

					case 2:	 iNeighbor = i - 1; jNeighbor = j; kNeighbor = k;
						temp = mesh_[level][i - 1][j][k];
						temp2 = mesh_[level][i - 2][j][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						break;
					case 10: iNeighbor = i - 1; jNeighbor = j;		kNeighbor = k + 1; break;
					case 11: iNeighbor = i - 1; jNeighbor = j;    	kNeighbor = k - 1; break;
					case 14: iNeighbor = i - 1; jNeighbor = j + 1;  kNeighbor = k;     break;
					case 17: iNeighbor = i - 1; jNeighbor = j - 1;  kNeighbor = k;     break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch

					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];

					//!check if neighbor was found
					if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
					mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
				} //!end dir
			} //!end k
		} //!end j

		//!do boundary nodes on inlet edges
		//!bottom inlet edge
		i = 0;
		int k = 0;

		for (int j = 0; j <= jmax; j++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				if (j == 0)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 0:
						//temp = mesh_[level][i + 1][j + 1][k + 1];
						//temp2 = mesh_[level][i + 2][j + 2][k + 2];
						//temp = mesh_[level][i + 1][j + 1][k];
						//temp2 = mesh_[level][i + 2][j + 2][k];
						temp = mesh_[level][i + 1][j][k];
						temp2 = mesh_[level][i + 2][j][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				else if (j == jmax)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 0:
						//temp = mesh_[level][i + 1][j - 1][k + 1];
						//temp2 = mesh_[level][i + 2][j - 2][k + 2];
						//temp = mesh_[level][i + 1][j - 1][k];
						//temp2 = mesh_[level][i + 2][j - 2][k];
						temp = mesh_[level][i + 1][j][k];
						temp2 = mesh_[level][i + 2][j][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				else
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 0:
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					case 6: iNeighbor = i + 1; jNeighbor = j; kNeighbor = k + 1;
						temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						temp2 = mesh_[level][iNeighbor + 1][jNeighbor][kNeighbor + 1];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						//!set neighbor
						neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end j

		//!top inlet edge
		k = kmax;
		for (int j = 0; j <= jmax; j++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				if (j == 0)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 0:
						//temp = mesh_[level][i + 1][j + 1][k - 1];
						//temp2 = mesh_[level][i + 2][j + 2][k - 2];
						//temp = mesh_[level][i + 1][j + 1][k];
						//temp2 = mesh_[level][i + 2][j + 2][k];
						temp = mesh_[level][i + 1][j][k];
						temp2 = mesh_[level][i + 2][j][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				if (j == jmax)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 0:
						//temp = mesh_[level][i + 1][j - 1][k - 1];
						//temp2 = mesh_[level][i + 2][j - 2][k - 2];
						//temp = mesh_[level][i + 1][j - 1][k];
						//temp2 = mesh_[level][i + 2][j - 2][k];
						temp = mesh_[level][i + 1][j][k];
						temp2 = mesh_[level][i + 2][j][k];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				else
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 0:
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					case 7: iNeighbor = i + 1; jNeighbor = j; kNeighbor = k - 1;
						temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						temp2 = mesh_[level][iNeighbor + 1][jNeighbor][kNeighbor - 1];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						//!set neighbor
						neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end j

		//!left inlet edge
		int j = 0;
		for (int k = 1; k <= kmax - 1; k++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
				{
				case 0:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 16:  iNeighbor = i + 1; jNeighbor = j + 1; kNeighbor = k;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor + 1][jNeighbor + 1][kNeighbor];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end k

		//!right inlet edge
		j = jmax;
		for (int k = 1; k <= kmax - 1; k++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
				{
				case 0:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 15:  iNeighbor = i + 1; jNeighbor = j - 1; kNeighbor = k;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor + 1][jNeighbor - 1][kNeighbor];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end k

		//!do boundary nodes on outlet edge
		//!bottom outlet edge
		i = imax;
		k = 0;
		for (int j = 0; j <= jmax; j++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				if (j == 0)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 4:
						temp = mesh_[level][i - 1][j + 1][k + 1];
						temp2 = mesh_[level][i - 2][j + 2][k + 2];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				else if (j == jmax)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 4:
						temp = mesh_[level][i - 1][j - 1][k + 1];
						temp2 = mesh_[level][i - 2][j - 2][k + 2];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				else
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 4:
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					case 10: iNeighbor = i - 1; jNeighbor = j; kNeighbor = k + 1;
						temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						temp2 = mesh_[level][iNeighbor - 1][jNeighbor][kNeighbor + 1];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						//!set neighbor
						neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end j

		//!top outlet edge
		k = kmax;
		for (int j = 0; j <= jmax; j++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				if (j == 0)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 5:
						temp = mesh_[level][i - 1][j + 1][k - 1];
						temp2 = mesh_[level][i - 2][j + 2][k - 2];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				if (j == jmax)
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 5:
						temp = mesh_[level][i - 1][j - 1][k - 1];
						temp2 = mesh_[level][i - 2][j - 2][k - 2];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}
				else
				{
					switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
					{
					case 5:
						mesh_[level][i][j][k]->setNormalDir(dir);
						continue;

					case 11: iNeighbor = i - 1; jNeighbor = j; kNeighbor = k - 1;
						temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						temp2 = mesh_[level][iNeighbor - 1][jNeighbor][kNeighbor - 1];
						mesh_[level][i][j][k]->setNormalNeighbor(temp);
						mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
						//!set neighbor
						neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
						break;

					default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
					} //!end switch
				}

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end j

		//!left outlet edge
		j = 0;
		for (int k = 1; k <= kmax - 1; k++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
				{
				case 1:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 14: iNeighbor = i - 1; jNeighbor = j + 1;	kNeighbor = k;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor - 1][jNeighbor + 1][kNeighbor];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end k

		//!right outlet edge
		j = jmax;

		for (int k = 1; k <= kmax - 1; k++)
		{
			//!corner/edge node?
			mesh_[level][i][j][k]->setIsCorner(true);

			for (int dir = 0; dir < 18; dir++)
			{
				Node* neighbor = NULL;
				Node* temp = NULL;
				Node* temp2 = NULL;
				int iNeighbor, jNeighbor, kNeighbor; //!node indices of next neighbor in this direction
				switch (dir) //!calc indices of neighbor --> set only pointers to existing non wall nodes; @bottom only nodes in top directions matter
				{
				case 3:
					mesh_[level][i][j][k]->setNormalDir(dir);
					continue;

				case 17: iNeighbor = i - 1; jNeighbor = j - 1; kNeighbor = k;
					temp = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					temp2 = mesh_[level][iNeighbor - 1][jNeighbor - 1][kNeighbor];
					mesh_[level][i][j][k]->setNormalNeighbor(temp);
					mesh_[level][i][j][k]->setNormalNeighborSnd(temp2);
					//!set neighbor
					neighbor = mesh_[level][iNeighbor][jNeighbor][kNeighbor];
					break;

				default: mesh_[level][i][j][k]->setNeighbor(dir, neighbor); continue; //!for other directions set neighbor to NULL and continue
				} //!end switch

				//!check if neighbor was found
				if (neighbor == NULL) cerr << " Error: could not find neighbor at node i j k" << i << " " << j << " " << k << ". Direction is " << dir;
				mesh_[level][i][j][k]->setNeighbor(dir, neighbor);	//!set pointer to neighbor
			} //!end dir
		} //!end k
	} //!end level
} //!end set connectivity


void Lattice::initialization(Control& bc) //!initialize all distribution values according to specified density and zero velocity
{
	int maxLevel = bc_.getLevelMax(); //!maximum level

	//!initialize populations based on zero-velocity equilibrium distribution
	for (int level = 0; level <= maxLevel; level++)		//!loop over all levels
	{
		double colFreq = (pow(bc_.getSpeedOfSound(), 2.) * (bc_.getTimeStep() / pow(2., level))) / (bc_.getKinematicViscosity() + 0.5 * pow(bc_.getSpeedOfSound(), 2.) * (bc_.getTimeStep() / pow(2., level))); //!calc dimensionless collision frequency
		double equil[19] = { 0. };
		double dens = bc_.getDensity();
		double timeStep = bc_.getTimeStep() / pow(2., level);
		double velo[3] = { 0. };
		bool   alternating = true;

		int imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
		int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
		int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level

		double Rc = 0.06;
		double A = 0.01;
		double xCenter = -8. * Rc;
		double yCenter = 0.;
		double cs = bc_.getSpeedOfSound();

		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction
				{
					mesh_[level][i][j][k]->setTimeStep(timeStep);

					//!initialize convected acoustic pulse
					double xNode = mesh_[level][i][j][k]->getX() - bc_.getChannelSizeX_() / 2.;
					double yNode = mesh_[level][i][j][k]->getY() - bc_.getChannelSizeY_() / 2.;
					double xRel = xNode - xCenter;
					double rhoNode = dens * (1 + A * exp(-(xRel*xRel + yNode*yNode) / (2 * Rc * Rc)));

                    if (bc.getOrderOfEquilibrium() == 2)
                    {
                        mesh_[level][i][j][k]->calcEquilibrium2_GHbasis(equil, bc, rhoNode, velo);
                    }
                    else if (bc.getOrderOfEquilibrium() == 3)
                    {
                        mesh_[level][i][j][k]->calcEquilibrium3_GHbasis(equil, bc, rhoNode, velo);
                    }
                    else if (bc.getOrderOfEquilibrium() == 4)
                    {
                        mesh_[level][i][j][k]->calcEquilibrium4_GHbasis(equil, bc, rhoNode, velo);
                    }
                    else
                    {
                        cerr << "Wrong order of equilibrium!" << flush;
                        exit(1);
                    }

					for (int dir = 0; dir <= 18; dir++) //!init all directions
					{
						mesh_[level][i][j][k]->setDistribution(true, dir, equil[dir]); //!set both distribution arrays
						mesh_[level][i][j][k]->setDistribution(false, dir, equil[dir]);
						mesh_[level][i][j][k]->setDistributionPreCol(dir, equil[dir]);
						mesh_[level][i][j][k]->setDistributionFictitious(dir, equil[dir]);
						mesh_[level][i][j][k]->setCubicMachCorrection(dir, 0.);
						mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, 0.);
						mesh_[level][i][j][k]->setCubicMachCorrectionHalf(dir, 0.);
						mesh_[level][i][j][k]->setCubicMachCorrectionFull(dir, 0.);
                        mesh_[level][i][j][k]->setCollisionFrequency(colFreq); //!set collision frequency

						double densityDummy = mesh_[level][i][j][k]->calcDensity(bc_, alternating);
						double veloDummy[3] = {0.};
						mesh_[level][i][j][k]->calcVelocity(bc_, alternating, densityDummy, veloDummy);

						mesh_[level][i][j][k]->setNodeDensity(densityDummy);
						int dirX = 0;
						int dirY = 1;
						int dirZ = 2;
						mesh_[level][i][j][k]->setNodeVelocityComponent(dirX, veloDummy[0]);
						mesh_[level][i][j][k]->setNodeVelocityComponent(dirY, veloDummy[1]);
						mesh_[level][i][j][k]->setNodeVelocityComponent(dirZ, veloDummy[2]);
					} //!end dir

					//!initialize strain rate and deviatoric stress tensor
					for (int x = 0; x < 3; x++)
					{
						for (int y = 0; y < 3; y++)
						{
							mesh_[level][i][j][k]->setStrainRateTensorFD(x, y, 0.);
							mesh_[level][i][j][k]->setStrainRateTensorFDfull(x, y, 0.);
							mesh_[level][i][j][k]->setStrainRateTensorFDhalf(x, y, 0.);
							mesh_[level][i][j][k]->setPi1Tensor(x, y, 0.);
						}
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end level
} //!end initialize


void Lattice::timeSteppingSRT(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

	//!safe pre-collision distribution on all nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				for (int k = 0; k <= kmax; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!do collision and streaming for all fluid nodes
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in y direction execpt walls
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction execpt walls
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternating); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternating); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternating); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternating);
				mesh_[level][i][jmax][k]->transport(negAlternating);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternating);
				mesh_[level][imax][j][k]->transport(negAlternating);
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				//!perform pressure outlet
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end j
	} //!end omp

	//!halfway bounce-back of all wall nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!perform bounce-back
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternating);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end i
	} //!end omp
} //!end timeSteppingSRT


void Lattice::timeSteppingMRT_GSbasis(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

	//!safe pre-collision distribution on all nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				for (int k = 0; k <= kmax; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!do collision and streaming for all fluid nodes
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in y direction execpt walls
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction execpt walls
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternating); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternating); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternating); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternating);
				mesh_[level][i][jmax][k]->transport(negAlternating);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternating);
				mesh_[level][imax][j][k]->transport(negAlternating);
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				//!perform pressure outlet
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end j
	} //!end omp

	//!halfway bounce-back of all wall nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!perform bounce-back
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternating);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end i
	} //!end omp
} //!end timeSteppingMRT_GSbasis


void Lattice::timeSteppingMRT_RMbasis(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution

	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

	//!safe pre-collision distribution on all nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				for (int k = 0; k <= kmax; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!do collision and streaming for all fluid nodes
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in y direction execpt walls
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction execpt walls
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternating); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternating); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternating); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternating);
				mesh_[level][i][jmax][k]->transport(negAlternating);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternating);
				mesh_[level][imax][j][k]->transport(negAlternating);
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				//!perform pressure outlet
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end j
	} //!end omp

	//!halfway bounce-back of all wall nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!perform bounce-back
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternating);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end i
	} //!end omp
} //!end timeSteppingMRT_RMbasis


void Lattice::timeSteppingRR(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution

	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

	//!safe pre-collision distribution on all nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				for (int k = 0; k <= kmax; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!do collision and streaming for all fluid nodes
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in y direction execpt walls
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction execpt walls
				{
					mesh_[level][i][j][k]->collideRR(bc_, alternating); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternating); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternating); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternating);
				mesh_[level][i][jmax][k]->transport(negAlternating);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternating);
				mesh_[level][imax][j][k]->transport(negAlternating);
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				//!perform pressure outlet
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end j
	} //!end omp

	//!halfway bounce-back of all wall nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!perform bounce-back
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternating);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end i
	} //!end omp
} //!end timeSteppingRR


void Lattice::timeSteppingHRR(int& level, bool& alternating) //!perform a collide stream step for this level
{
	//!the boolean variable alternating is used to differentiate between the two distribution arrays (0/1)
	bool negAlternating = !alternating; //!accordingly by negAlternating is opposite array of distribution
	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

	//!safe pre-collision distribution on all nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				for (int k = 0; k <= kmax; k++) //!nodes in z direction from coarse lower interface to upper interface
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternating, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collide and stream
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!do collision and streaming for all fluid nodes
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in y direction execpt walls
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction execpt walls
				{
					mesh_[level][i][j][k]->collideHRR(bc_, alternating); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternating); //!transport the distribution from alternating (just performed collide) to neighbors (write into negAlternating array)
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternating); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternating); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternating);
				mesh_[level][i][jmax][k]->transport(negAlternating);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternating);
				mesh_[level][imax][j][k]->transport(negAlternating);
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				//!perform pressure outlet
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end j
	} //!end omp

	//!halfway bounce-back of all wall nodes
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!perform bounce-back
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternating);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternating);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternating);
			} //!end k
		} //!end i
	} //!end omp
} //!end timeSteppingHRR


void Lattice::solve(bool& alternatingL0)
{
	int level = 0;
	string collisionModel = bc_.getCollisionModel();

	if (collisionModel == "SRT")
	{
		this->timeSteppingSRT(level, alternatingL0);
	}
	else if (collisionModel == "MRT-GS")
	{
		this->timeSteppingMRT_GSbasis(level, alternatingL0);
	}
	else if (collisionModel == "MRT-RM")
	{
		this->timeSteppingMRT_RMbasis(level, alternatingL0);
	}
	else if (collisionModel == "RR")
	{
		this->timeSteppingRR(level, alternatingL0);
	}
	else if (collisionModel == "HRR")
	{
		this->timeSteppingHRR(level, alternatingL0);
	}
} //!end solve

void Lattice::solve(bool& alternatingL0, bool& alternatingL1)
{
	string collisionModel = bc_.getCollisionModel();

	if (collisionModel == "SRT")
	{
		this->nestedTimeSteppingSRT(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "MRT-GS")
	{
		this->nestedTimeSteppingMRT_GSbasis(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "MRT-RM")
	{
		this->nestedTimeSteppingMRT_RMbasis(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "RR")
	{
		this->nestedTimeSteppingRR(alternatingL1, alternatingL0, bc_);
	}
	else if (collisionModel == "HRR")
	{
		this->nestedTimeSteppingHRR(alternatingL1, alternatingL0, bc_);
	}
} //!end solve


void Lattice::writeResultsVTK(string filename, bool& alternatingL0, bool& alternatingL1)
{
	filename.append(".vtk");
	ofstream paraview(filename.c_str()); //!open file stream to write

	if (!paraview) //!if stream is corrupted, quit
	{
		cerr << "An error occured while writing paraview file\n";
		return;
	}
	else
	{
		//!write header of paraview file
		paraview << "# vtk DataFile Version 2.0" << "\n";
		paraview << "LBM results for grid refinement benchmark in channel flow\n";
		paraview << "ASCII" << "\n";
		paraview << "DATASET UNSTRUCTURED_GRID" << "\n";

		paraview << "POINTS ";
		int pntPos = static_cast<int>(paraview.tellp()); //!remember stream position to write number of nodes to the file (this will get important later)
		paraview << "              " << " float" << "\n";
		paraview << setiosflags(ios::left | ios::showpoint | ios::scientific) << setprecision(8);

		int nodeCntrBot = 0; //!number of fine nodes in bottom region
		int nodeCntrTop = 0; //!number of fine nodes in top region
		int nodeCntrLft = 0; //!number of fine nodes in left region
		int nodeCntrRgt = 0; //!number of fine nodes in right region
		int nodeCntrFnt = 0; //!number of fine nodes in inlet region
		int nodeCntrBck = 0; //!number of fine nodes in outlet region
		int nodeCntrMid = 0; //!number of coarse nodes in mid region
		int pntCntr = 0;

		if (bc_.getLevelMax() > 0)
		{
			//------------------------------------------//
			//        coordinates in *.vtk-file:		//
			//------------------------------------------//

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			int level = 1;
			int imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			int jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			int kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level

			int refLayer = bc_.getRefinementLayer() * pow(2, level);
			int refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrBot++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrTop++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction
			{
				for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrLft++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction except those considered before
			{
				for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrRgt++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction except those considered before
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction except those considered before
				{
					for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrFnt++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction except those considered before
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction except those considered before
				{
					for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrBck++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX();

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in unrefined core region
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in unrefined core region
				{
					for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction in unrefined core region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							nodeCntrMid++; //!count nodes
						}
					} //!end i
				} //!end j
			} //!end k

			//!assemble and assign index
			paraview << "\n";
			int nodeCntrTot = nodeCntrBot + nodeCntrTop + nodeCntrLft + nodeCntrRgt + nodeCntrFnt + nodeCntrBck; //!total number of fine nodes
			pntCntr = nodeCntrTot + nodeCntrMid; //!store total number of nodes

			paraview << "CELLS "; //!write cells
			int cellPos = static_cast<int>(paraview.tellp()); //!remember cell position
			int cellCntr = 0;
			paraview << "         " << " " << "         " << "\n";

			//------------------------------------------//
			//		  connectivity in *.vtk-file		//
			//------------------------------------------//

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			int icmax = imax - 1; //!index of last cell in x direction in ->this level
			int jcmax = jmax - 1; //!index of last cell in y direction in ->this level
			int kcmax = kmax - 1; //!index of last cell in z direction in ->this level

			int nodeCntrX = imax + 1; //!number of nodes in x direction in this level
			int nodeCntrY = jmax + 1; //!number of nodes in y direction in this level
			int nodeCntrXY = nodeCntrX * nodeCntrY; //!number of nodes per xy plane in this level
			int nodeCntrRef = refLayer + 1;
			int nodeCntrRefX = refLayerX + 1;
			int nodeCntrZ = 0;
			int nodeCntrJ = 0;
			int nodeCntrI = 0;

			for (int k = 0; k <= refLayer - 1; k++) //!all cells in z direction in refined bottom region
			{
				for (int j = 0; j <= jcmax; j++) //!all cells in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + j * nodeCntrX + k * nodeCntrXY;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrXY << " " << firstNodeIndex + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrXY + nodeCntrX << " " << firstNodeIndex + nodeCntrXY + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kcmax - refLayer + 1; k <= kcmax; k++) //!all cells in z direction in refined top region
			{
				for (int j = 0; j <= jcmax; j++) //!all cells in y direction
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + j * nodeCntrX + nodeCntrZ * nodeCntrXY + nodeCntrBot;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrXY << " " << firstNodeIndex + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrXY + nodeCntrX << " " << firstNodeIndex + nodeCntrXY + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells

						}
					} //!end i
				} //!end j
				nodeCntrZ++;
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			nodeCntrZ = 0;

			for (int k = refLayer; k <= kcmax - refLayer; k++) //!all cells in z direction in refined left region except cells already considered before
			{
				for (int j = 0; j <= refLayer - 1; j++) //!all cells in y direction in refined left region
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + j * nodeCntrX + nodeCntrZ * nodeCntrRef * nodeCntrX + nodeCntrBot + nodeCntrTop;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
				} //!end j
				nodeCntrZ++;
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			nodeCntrZ = 0;

			for (int k = refLayer; k <= kcmax - refLayer; k++) //!all cells in z direction in refined right region except nodes already considered before
			{
				for (int j = jcmax - refLayer + 1; j <= jcmax; j++) //!all cells in y direction in refined right region
				{
					for (int i = 0; i <= icmax; i++) //!all cells in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + nodeCntrJ * nodeCntrX + nodeCntrZ * nodeCntrRef * nodeCntrX + nodeCntrBot + nodeCntrTop + nodeCntrLft;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX << " " << firstNodeIndex + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX << " " << firstNodeIndex + nodeCntrRef * nodeCntrX + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
					nodeCntrJ++;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			nodeCntrZ = 0;
			nodeCntrJ = 0;

			for (int k = refLayer; k <= kcmax - refLayer; k++) //!all cells in z direction in refined inlet region except nodes already considered before
			{
				for (int j = refLayer; j <= jcmax - refLayer; j++) //!all cells in y direction in refined inlet region except nodes already considered before
				{
					for (int i = 0; i <= refLayerX - 1; i++) //!all cells in x direction in refined inlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = i + nodeCntrJ * nodeCntrRefX + nodeCntrZ * nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrBot + nodeCntrTop + nodeCntrLft + nodeCntrRgt;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrRefX << " " << firstNodeIndex + nodeCntrRefX + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) << " " << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRefX << " " << firstNodeIndex + nodeCntrRefX * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRefX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
					} //!end i
					nodeCntrJ++;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			nodeCntrZ = 0;
			nodeCntrJ = 0;

			for (int k = refLayer; k <= kcmax - refLayer; k++) //!all cells in z direction in refined outlet region except nodes already considered before
			{
				for (int j = refLayer; j <= jcmax - refLayer; j++) //!all cells in y direction in refined outlet region except nodes already considered before
				{
					for (int i = icmax - refLayer + 1; i <= icmax; i++) //!all cells in x direction in refined outlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = nodeCntrI + nodeCntrJ * nodeCntrRef + nodeCntrZ * nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrBot + nodeCntrTop + nodeCntrLft + nodeCntrRgt + nodeCntrFnt;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrRef << " " << firstNodeIndex + nodeCntrRef + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) << " " << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRef << " " << firstNodeIndex + nodeCntrRef * (nodeCntrY - 2 * (nodeCntrRef - 1)) + nodeCntrRef + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cells
						}
						nodeCntrI++;
					} //!end i
					nodeCntrJ++;
					nodeCntrI = 0;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k


			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX();

			icmax = imax - 1; //!index of last cell in x direction in ->this level
			jcmax = jmax - 1; //!index of last cell in y direction in ->this level
			kcmax = kmax - 1; //!index of last cell in z direction in ->this level

			nodeCntrRef = refLayer;
			nodeCntrRefX = refLayerX;

			nodeCntrX = imax + 1; //!number of nodes in x direction in ->this level
			nodeCntrY = jmax + 1; //!number of nodes in y direction in ->this level
			nodeCntrXY = (nodeCntrX - (nodeCntrRefX + nodeCntrRef)) * (nodeCntrY - 2 * nodeCntrRef); //!nodes per x y plane in unrefined core

			paraview << "\n";

			nodeCntrZ = 0;
			nodeCntrJ = 0;
			nodeCntrI = 0;

			for (int k = refLayer; k <= kcmax - refLayer; k++) //!all cells in z direction in unrefined mid region
			{
				for (int j = refLayer; j <= jcmax - refLayer; j++) //!all cells in y direction in unrefined mid region
				{
					for (int i = refLayerX; i <= icmax - refLayer; i++) //!all cells in x direction in unrefined mid region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int firstNodeIndex = nodeCntrI + nodeCntrJ * (nodeCntrX - (nodeCntrRefX + nodeCntrRef)) + nodeCntrZ * nodeCntrXY + nodeCntrTot;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << firstNodeIndex << " " << firstNodeIndex + 1 << " " << firstNodeIndex + nodeCntrX - (nodeCntrRefX + nodeCntrRef) << " " << firstNodeIndex + nodeCntrX - (nodeCntrRefX + nodeCntrRef) + 1 << " "; //!nodes 0 to 3
							paraview << firstNodeIndex + nodeCntrXY << " " << firstNodeIndex + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << firstNodeIndex + nodeCntrXY + nodeCntrX - (nodeCntrRefX + nodeCntrRef) << " " << firstNodeIndex + nodeCntrXY + nodeCntrX - (nodeCntrRefX + nodeCntrRef) + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cell
						}
						nodeCntrI++;
					} //!end i
					nodeCntrJ++;
					nodeCntrI = 0;
				} //!end j
				nodeCntrZ++;
				nodeCntrJ = 0;
			} //!end k

			paraview << "\n";
			paraview << "CELL_TYPES " << cellCntr << "\n"; //!write cell types to file

			for (int i = 0; i < cellCntr; i++)		//!write cell types to paraview file
			{
				paraview << "11" << "\n";
			}
			paraview << "\n";

			//-------------------------//
			// solutions in *.vtk-file //
			//-------------------------//

			//->> COLLISION FREQUENCY <<-//

			paraview << "POINT_DATA " << pntCntr << "\n";
			paraview << "SCALARS Omega float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined left region except nodes already considered before
			{
				for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined right region except nodes already considered before
			{
				for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined inlet region except nodes already considered before
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
				{
					for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
				{
					for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX();

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in unrefined mid region
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k


			//->> NODE ID <<-//

			paraview << "SCALARS NodeID float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getNodeID() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// top part level 1 //
			//------------------//

			for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getNodeID() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//-------------------//
			// left part level 1 //
			//-------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined left region except nodes already considered before
			{
				for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getNodeID() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// right part level 1 //
			//--------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined right region except nodes already considered before
			{
				for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getNodeID() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined inlet region except nodes already considered before
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
				{
					for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getNodeID() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
				{
					for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								paraview << mesh_[level][i][j][k]->getNodeID() << "\n";
							}
						}
					} //!end i
				} //!end j
			} //!end k

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX();

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in unrefined mid region
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getNodeID() << "\n";
						}
					} //!end i
				} //!end j
			} //!end k


			//->> PRESSURE <<-//

			paraview << "SCALARS Pressure float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";
			double speedOfSoundSq = pow(bc_.getSpeedOfSound(), 2.); //!square of speed of sound

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//------------------//
			// top part level 1 //
			//------------------//
			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//-------------------//
			// left part level 1 //
			//-------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//--------------------//
			// right part level 1 //
			//--------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) * speedOfSoundSq << "\n"; //!for fluid nodes calculate pressure
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes write pressure calculated from specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes write pressure calculated from specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) * speedOfSoundSq << "\n"; //!for boundary nodes @ inlet get pressure from +x neighbor
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
								{
									if (mesh_[level][i][j][k]->getTag() == "wall")
									{
										paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write pressure calculated from specified density
									}
									else if (mesh_[level][i][j][k]->getTag() == "boundary")
									{
										bool negAlternatingL1 = !alternatingL1;
										paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) * speedOfSoundSq << "\n";
									}
									else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX();

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in unrefined mid region
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
						}
					}//end i
				}//end j
			}//end k


			//->> Density <<-//

			paraview << "SCALARS Density float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) << "\n"; //!for boundary nodes @ inlet get density from +x neighbor
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//------------------//
			// top part level 1 //
			//------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) << "\n"; //!for boundary nodes @ inlet get density from +x neighbor
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//-------------------//
			// left part level 1 //
			//-------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) << "\n"; //!for boundary nodes @ inlet get density from +x neighbor
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//--------------------//
			// right part level 1 //
			//--------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) << "\n"; //!for boundary nodes @ inlet get density from +x neighbor
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes write specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes write specified density
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlterantingL1 = !alternatingL1;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlterantingL1) << "\n"; //!for boundary nodes @ inlet get density from +x neighbor
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
								{
									if (mesh_[level][i][j][k]->getTag() == "wall")
									{
										paraview << bc_.getDensity() << "\n"; //!for wall nodes or boundary nodes @ outlet write specified density
									}
									else if (mesh_[level][i][j][k]->getTag() == "boundary")
									{
										bool negAlternatingL1 = !alternatingL1;
										paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1) << "\n";
									}
									else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1) << "\n"; //!for fluid and boundary nodes calc density
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX();

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in unrefined mid region
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0) << "\n"; //!for fluid and boundary nodes calc density
						}
					}//end i
				}//end j
			}//end k

			//->> VELOCITY <<-//

			paraview << "\nVECTORS Velocity float " << "\n";

			//---------------------//
			// bottom part level 1 //
			//---------------------//

			level = 1;
			imax = (bc_.getNumberOfNodesX() - 1) * pow(2, level); //!index of last node in x direction in ->this level
			jmax = (bc_.getNumberOfNodesY() - 1) * pow(2, level); //!index of last node in y direction in ->this level
			kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2, level); //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer() * pow(2, level);
			refLayerX = bc_.getRefinementLayerX() * pow(2, level);

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								//else if ((mesh_[level][i][j][k]->getTag() == "boundary") && (i == imax))
								//{
								//	double dens = mesh_[level][i][j][k]->getNeighbor(2)->calcDensity(bc_, alternatingL1);
								//	mesh_[level][i][j][k]->getNeighbor(2)->calcVelocity(bc_, alternatingL1, dens, velo);
								//	paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for boundary nodes @ outlet get velocity from -x neighbor
								//}
								//else if ((mesh_[level][i][j][k]->getTag() == "boundary") && (i == 0))
								//{
								//	paraview << bc_.getVelocityInlet()[0] << " " << bc_.getVelocityInlet()[1] << " " << bc_.getVelocityInlet()[2] << " " << "\n"; //!for boundary nodes @ inlet write specified velocity
								//}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = 0; k <= refLayer; k++) //!all nodes in z direction in refined bottom region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, negAlternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for boundary nodes @ outlet get velocity from -x neighbor
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//------------------//
			// top part level 1 //
			//------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = kmax - refLayer; k <= kmax; k++) //!all nodes in z direction in refined top region
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, negAlternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for boundary nodes @ outlet get velocity from -x neighbor
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//-------------------//
			// left part level 1 //
			//-------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = 0; j <= refLayer; j++) //!all nodes in y direction in refined left region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, negAlternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for boundary nodes @ outlet get velocity from -x neighbor
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//--------------------//
			// right part level 1 //
			//--------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = jmax - refLayer; j <= jmax; j++) //!all nodes in y direction in refined right region
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, negAlternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for boundary nodes @ outlet get velocity from -x neighbor
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//--------------------//
			// inlet part level 1 //
			//--------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = 0; i <= refLayerX; i++) //!all nodes in x direction in refined inlet region
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall")
								{
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
								}
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL1 = !alternatingL1;
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, negAlternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
								}
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//---------------------//
			// outlet part level 1 //
			//---------------------//

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
						{
							double velo[3] = { 0. };
							if (mesh_[level][i][j][k]->getTag() == "wall")
							{
								paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
							}
							else
							{
								double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
								mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
								paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
							}
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in refined outlet region except nodes already considered before
				{
					for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction in refined inlet region except nodes already considered before
					{
						for (int i = imax - refLayer; i <= imax; i++) //!all nodes in x direction in refined outlet region
						{
							double velo[3] = { 0. };
							if (mesh_[level][i][j][k]->getTag() == "wall")
							{
								paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at walls
							}
							else if (mesh_[level][i][j][k]->getTag() == "boundary")
							{
								bool negAlternatingL1 = !alternatingL1;
								double dens = mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL1);
								mesh_[level][i][j][k]->calcVelocity(bc_, negAlternatingL1, dens, velo);
								paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for boundary nodes @ outlet get velocity from -x neighbor
							}
							else
							{
								double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL1);
								mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL1, dens, velo);
								paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!for fluid nodes calc velocity
							}
						} //!end i
					} //!end j
				} //!end k
			}

			//------------------//
			// mid part level 0 //
			//------------------//

			level = 0;
			imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
			refLayer = bc_.getRefinementLayer();
			refLayerX = bc_.getRefinementLayerX();

			for (int k = refLayer; k <= kmax - refLayer; k++) //!all nodes in z direction in unrefined mid region
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
					{
						double velo[3] = { 0. };
						double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0);
						mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL0, dens, velo);
						paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
					} //!end i
				} //!end j
			} //!end k

			//-------//
			// final //
			//-------//

			//!update node and cell numbers
			paraview.seekp(pntPos); //!go to nodePosition of the stream pointer
			paraview << pntCntr;	//!insert correct node number
			paraview.seekp(cellPos); //!go to cellPosition of the stream pointer
			paraview << cellCntr << " " << 9 * cellCntr; //!insert correct cell number and number of indices
			paraview.close(); //!close file
		} //!end multi-level

		//!SINGLE LEVEL
		else
		{
			//!write all node coordinates to file
			int level = 0;
			int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
			int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
			int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

			int nodeCntrX = imax + 1; //!nodes in x direction in this level
			int nodeCntrY = jmax + 1; //!nodes in y direction in this level
			int nodeCntrXY = nodeCntrX * nodeCntrY; //!nodes per xy plane in this level

			int index = 0;

			for (int k = 0; k <= kmax; k++) //!all nodes in z direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getX() << " " << mesh_[level][i][j][k]->getY() << " " << mesh_[level][i][j][k]->getZ() << "\n"; //!write coordinates
							index++; //!count nodes
						} //!end if
					} //!end i
				} //!end j
			} //!end k

			paraview << "\n";

			pntCntr = index;	//!store numberOfPoints

			paraview << "CELLS "; //!write cells
			int cellPos = static_cast<int>(paraview.tellp()); //!remember cell position
			int cellCntr = 0;
			paraview << "         " << " " << "         " << "\n";

			//!write cell connectivity and write file; needed for visualization

			for (int k = 0; k < kmax; k++) //!all nodes in z direction
			{
				for (int j = 0; j < jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i < imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!cells can only exist if node exists
						{
							//!write indices of nodes, who create the cell to file
							//!paraview cell connectivity (for type vtk_voxel (=11)) is:
							/*
							  6--------7
							 /|	      /|
							4-|------5 |
							| 2------|-3   z  y
							|/		 |/	    |/
							0--------1		-->x
							*/
							//!calc index of node 0 (i.e. line in paraview file) from i,j and k
							int nodeIndex0 = i + j * nodeCntrX + k * nodeCntrXY;

							//!write connectivity
							paraview << "8 "; //!write number of points in cell (i.e. number of variables to read in this line)
							paraview << nodeIndex0 << " " << nodeIndex0 + 1 << " " << nodeIndex0 + nodeCntrX << " " << nodeIndex0 + nodeCntrX + 1 << " "; //!nodes 0 to 3
							paraview << nodeIndex0 + nodeCntrXY << " " << nodeIndex0 + nodeCntrXY + 1 << " "; //!nodes 4 and 5
							paraview << nodeIndex0 + nodeCntrXY + nodeCntrX << " " << nodeIndex0 + nodeCntrXY + nodeCntrX + 1 << "\n"; //!nodes 6 and 7
							cellCntr++; //!count cell
						} //!end if
					} //!end i
				} //!end j
			} //!end k

			paraview << "\n";
			paraview << "CELL_TYPES " << cellCntr << "\n"; //!write cell types to file

			for (int i = 0; i < cellCntr; i++)		//!write cell types to paraview file
			{
				paraview << "11" << "\n";
			}
			paraview << "\n";

			//!write all node data to file

			//!first of all: pressure
			paraview << "POINT_DATA " << pntCntr << "\n";
			paraview << "SCALARS Pressure float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			double speedOfSoundSq = pow(bc_.getSpeedOfSound(), 2.); //!square of speed of sound

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall") paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes write pressure calculated from specified density
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							} //!end if
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall") paraview << speedOfSoundSq * bc_.getDensity() << "\n"; //!for wall nodes write pressure calculated from specified density
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL0 = !alternatingL0;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL0) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0) * speedOfSoundSq << "\n"; //!for fluid and boundary nodes calc pressure
							} //!end if
						} //!end i
					} //!end j
				} //!end k
			}

			//!density
			paraview << "SCALARS Density float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall") paraview << bc_.getDensity() << "\n"; //!for wall nodes write specified density
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0) << "\n"; //!for fluid and boundary nodes calc density
							} //!end if
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								if (mesh_[level][i][j][k]->getTag() == "wall") paraview << bc_.getDensity() << "\n"; //!for wall nodes write specified density
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL0 = !alternatingL0;
									paraview << mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL0) << "\n"; //!for fluid and boundary nodes calc density
								}
								else paraview << mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0) << "\n"; //!for fluid and boundary nodes calc density
							} //!end if
						} //!end i
					} //!end j
				} //!end k
			}

			paraview << "SCALARS Omega float 1" << "\n";
			paraview << "LOOKUP_TABLE default" << "\n";

			for (int k = 0; k <= kmax; k++) //!all nodes in z direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					for (int i = 0; i <= imax; i++) //!all nodes in x direction
					{
						if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
						{
							paraview << mesh_[level][i][j][k]->getCollisionFrequency() << "\n";
						} //!end if
					} //!end i
				} //!end j
			} //!end k

			//!next: write velocity to file
			paraview << "\nVECTORS Velocity float " << "\n";

			if (bc_.getReflectionHandling() == "no")
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall") paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at wall nodes
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL0, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
								}
							} //!end if
						} //!end i
					} //!end j
				} //!end k
			}
			else
			{
				for (int k = 0; k <= kmax; k++) //!all nodes in z direction
				{
					for (int j = 0; j <= jmax; j++) //!all nodes in y direction
					{
						for (int i = 0; i <= imax; i++) //!all nodes in x direction
						{
							if (mesh_[level][i][j][k] != NULL) //!write only existing nodes
							{
								double velo[3] = { 0. };
								if (mesh_[level][i][j][k]->getTag() == "wall") paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!no slip condition at wall nodes
								else if ((i > 0) && (mesh_[level][i][j][k]->getTag() == "boundary"))
								{
									bool negAlternatingL0 = !alternatingL0;
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, negAlternatingL0);
									mesh_[level][i][j][k]->calcVelocity(bc_, negAlternatingL0, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
								}
								else
								{
									double dens = mesh_[level][i][j][k]->calcDensity(bc_, alternatingL0);
									mesh_[level][i][j][k]->calcVelocity(bc_, alternatingL0, dens, velo);
									paraview << velo[0] << " " << velo[1] << " " << velo[2] << " " << "\n"; //!write velocity to file
								}
							} //!end if
						} //!end i
					} //!end j
				} //!end k
			}

			//!update node and cell numbers
			paraview.seekp(pntPos);	//!go to nodePosition of the stream pointer
			paraview << pntCntr;	//!insert correct node number
			paraview.seekp(cellPos); //!go to cellPosition of the stream pointer
			paraview << cellCntr << " " << 9 * cellCntr;	//!insert correct cell number and number of indices
			paraview.close(); //!close file
		} //!end else
	} //!end else
} //!end write results


void Lattice::nestedTimeSteppingSRT(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	int d1Interpol; //!switch for 1-D interpolation in different cartesian directions
	int d2Interpol; //!switch for 2-D interpolation in different cartesian directions
	int refLayer; //!number of refined cells
	int refLayerX; //!number of refined cells
	int overlap; //!number of coarse overlap-cells

	//!boolean switches to change distribution arrays
	bool negAlternatingL1 = !alternatingL1;
	bool negAlternatingL0 = !alternatingL0;

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             collision and transport				//
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!current distribution (t=0) level 0: alternating
	//!current distribution (t=0) level 1: alternating

	//!level 0
	int level = 0;
	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined cells
	refLayerX = bc_.getRefinementLayerX(); //!number of refined cells
	overlap = bc_.getOverlap(); //!number of coarse overlap-cells

	/*SAFE PRE-COLLISION DISTRIBUTIONS ON COARSE NODES FOR TEMPORAL INTERPOLATION*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collision and transport on level 0
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternatingL0); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL0);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!level 1
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
	jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
	kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

	//!collision and transport on level 1 (relevant nodes only)
	//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!bottom refined region
				for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!top refined region
				for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!left refined region
				for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!right refined region
				for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!left refined region
				for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!right refined region
				for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternatingL1); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternatingL1); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternatingL1); //!left bc
				mesh_[level][i][jmax][k]->transport(negAlternatingL1); //!right bc
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternatingL1);
				mesh_[level][imax][j][k]->transport(negAlternatingL1);
			} //!end k
		} //!end j
	} //!end omp

	//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end k
		} //!end j
	} //!end omp

		//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!bottom and top wall
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	// temporal, spatial interpolation and coupling c2f //
	//--------------------------------------------------//
	//--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
			{
				//!bottom interface
				int k = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!top interface
				k = kmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!right interface
				j = jmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!inlet interface
				int i = refLayerX;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!outlet interface
				i = imax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for z = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x y plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				//!interpolation in x y plane
				d2Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for y = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in x z plane
				d2Interpol = 2;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for x = const. //
	//--------------------------------------//

	//--------------------------//
	// mid-nodes in y direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------//
	// mid-nodes in z direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//---------------------------//
	// center-nodes in y z plane //
	//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in y z plane
				d2Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	//        collision and transport (level 1)			//
	//--------------------------------------------------//
	//--------------------------------------------------//

//!collision and transport on level 1 (relevant nodes only)
//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!bottom refined region
				for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!top refined region
				for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!left refined region
				for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!right refined region
				for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!left refined region
				for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!right refined region
				for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
				{
					mesh_[level][i][j][k]->collideSRT(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(alternatingL1); //!bottom bc
				mesh_[level][i][j][kmax]->transport(alternatingL1); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(alternatingL1); //!left bc
				mesh_[level][i][jmax][k]->transport(alternatingL1); //!right bc
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(alternatingL1);
				mesh_[level][imax][j][k]->transport(alternatingL1);
			} //!end k
		} //!end j
	} //!end omp

	//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end k
		} //!end j
	} //!end omp

		//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!bottom and top wall
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             fine-to-coarse coupling			    //
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!level 0
	level = 0;
	refLayer = bc_.getRefinementLayer();
	refLayerX = bc_.getRefinementLayerX();
	imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				//!bottom interface
				int k = refLayer - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}

				//!top interface
				k = kmax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}

				//!right interface
				j = jmax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!y direction
		{
			for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!z direction
			{
				//!inlet interface
				int i = refLayerX - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}

				//!outlet interface
				i = imax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
			} //!end k
		} //!end j
	} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             coarse-to-fine coupling			     //
		//--------------------------------------------------//
		//--------------------------------------------------//

	//!level 1
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
	jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
	kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
			{
				//!bottom interface
				int k = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!top interface
				k = kmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!right interface
				j = jmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!inlet interface
				int i = refLayerX;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!outlet interface
				i = imax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for z = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x y plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				//!interpolation in x y plane
				d2Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for y = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in x z plane
				d2Interpol = 2;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for x = const. //
	//--------------------------------------//

	//--------------------------//
	// mid-nodes in y direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------//
	// mid-nodes in z direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//---------------------------//
	// center-nodes in y z plane //
	//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in y z plane
				d2Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp
} //!end nestedTimeSteppingSRT


void Lattice::nestedTimeSteppingMRT_GSbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	int d1Interpol; //!switch for 1-D interpolation in different cartesian directions
	int d2Interpol; //!switch for 2-D interpolation in different cartesian directions
	int refLayer; //!number of refined cells
	int refLayerX; //!number of refined cells
	int overlap; //!number of coarse overlap-cells

	//!boolean switches to change distribution arrays
	bool negAlternatingL1 = !alternatingL1;
	bool negAlternatingL0 = !alternatingL0;

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             collision and transport				//
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!current distribution (t=0) level 0: alternating
	//!current distribution (t=0) level 1: alternating

	//!level 0
	int level = 0;
	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined cells
	refLayerX = bc_.getRefinementLayerX(); //!number of refined cells
	overlap = bc_.getOverlap(); //!number of coarse overlap-cells

	/*SAFE PRE-COLLISION DISTRIBUTIONS ON COARSE NODES FOR TEMPORAL INTERPOLATION*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collision and transport on level 0
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL0); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL0);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!level 1
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
	jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
	kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

	//!collision and transport on level 1 (relevant nodes only)
	//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!bottom refined region
				for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!top refined region
				for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!left refined region
				for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!right refined region
				for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!left refined region
				for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!right refined region
				for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternatingL1); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternatingL1); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternatingL1); //!left bc
				mesh_[level][i][jmax][k]->transport(negAlternatingL1); //!right bc
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternatingL1);
				mesh_[level][imax][j][k]->transport(negAlternatingL1);
			} //!end k
		} //!end j
	} //!end omp

	//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end k
		} //!end j
	} //!end omp

		//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!bottom and top wall
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	// temporal, spatial interpolation and coupling c2f //
	//--------------------------------------------------//
	//--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
			{
				//!bottom interface
				int k = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!top interface
				k = kmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!right interface
				j = jmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!inlet interface
				int i = refLayerX;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!outlet interface
				i = imax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for z = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x y plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				//!interpolation in x y plane
				d2Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for y = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in x z plane
				d2Interpol = 2;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for x = const. //
	//--------------------------------------//

	//--------------------------//
	// mid-nodes in y direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------//
	// mid-nodes in z direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//---------------------------//
	// center-nodes in y z plane //
	//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in y z plane
				d2Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	//        collision and transport (level 1)			//
	//--------------------------------------------------//
	//--------------------------------------------------//

//!collision and transport on level 1 (relevant nodes only)
//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!bottom refined region
				for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!top refined region
				for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!left refined region
				for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!right refined region
				for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!left refined region
				for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!right refined region
				for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
				{
					mesh_[level][i][j][k]->collideMRT_GSbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(alternatingL1); //!bottom bc
				mesh_[level][i][j][kmax]->transport(alternatingL1); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(alternatingL1); //!left bc
				mesh_[level][i][jmax][k]->transport(alternatingL1); //!right bc
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(alternatingL1);
				mesh_[level][imax][j][k]->transport(alternatingL1);
			} //!end k
		} //!end j
	} //!end omp

	//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end k
		} //!end j
	} //!end omp

		//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!bottom and top wall
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             fine-to-coarse coupling			    //
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!level 0
	level = 0;
	refLayer = bc_.getRefinementLayer();
	refLayerX = bc_.getRefinementLayerX();
	imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				//!bottom interface
				int k = refLayer - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}

				//!top interface
				k = kmax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}

				//!right interface
				j = jmax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!y direction
		{
			for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!z direction
			{
				//!inlet interface
				int i = refLayerX - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}

				//!outlet interface
				i = imax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
					}
				}
			} //!end k
		} //!end j
	} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             coarse-to-fine coupling			     //
		//--------------------------------------------------//
		//--------------------------------------------------//

	//!level 1
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
	jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
	kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
			{
				//!bottom interface
				int k = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!top interface
				k = kmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!right interface
				j = jmax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!inlet interface
				int i = refLayerX;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

				//!outlet interface
				i = imax - refLayer;
				if (bc_.getScalingStyle() == "SRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getScalingStyle() == "MRTSTYLE")
				{
					mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
				}

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for z = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x y plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				//!interpolation in x y plane
				d2Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for y = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in x z plane
				d2Interpol = 2;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for x = const. //
	//--------------------------------------//

	//--------------------------//
	// mid-nodes in y direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------//
	// mid-nodes in z direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//---------------------------//
	// center-nodes in y z plane //
	//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in y z plane
				d2Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp
} //!end nestedTimeSteppingMRT_GSbasis


void Lattice::nestedTimeSteppingMRT_RMbasis(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	int d1Interpol; //!switch for 1-D interpolation in different cartesian directions
	int d2Interpol; //!switch for 2-D interpolation in different cartesian directions
	int refLayer; //!number of refined cells
	int refLayerX; //!number of refined cells
	int overlap; //!number of coarse overlap-cells

	//!boolean switches to change distribution arrays
	bool negAlternatingL1 = !alternatingL1;
	bool negAlternatingL0 = !alternatingL0;

	//--------------------------------------------------//
	//--------------------------------------------------//
	//             collision and transport				//
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!current distribution (t=0) level 0: alternating
	//!current distribution (t=0) level 1: alternating

	//!level 0
	int level = 0;
	int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
	refLayer = bc_.getRefinementLayer(); //!number of refined cells
	refLayerX = bc_.getRefinementLayerX(); //!number of refined cells
	overlap = bc_.getOverlap(); //!number of coarse overlap-cells

	/*SAFE PRE-COLLISION DISTRIBUTIONS ON COARSE NODES FOR TEMPORAL INTERPOLATION*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
				{
					for (int dir = 0; dir < 19; dir++)
					{
						double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
						mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
					}
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!collision and transport on level 0
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL0); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL0);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!level 1
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
	jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
	kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

	//!collision and transport on level 1 (relevant nodes only)
	//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!bottom refined region
				for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!top refined region
				for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!left refined region
				for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!right refined region
				for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!left refined region
				for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!right refined region
				for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, alternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(negAlternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(negAlternatingL1); //!bottom bc
				mesh_[level][i][j][kmax]->transport(negAlternatingL1); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(negAlternatingL1); //!left bc
				mesh_[level][i][jmax][k]->transport(negAlternatingL1); //!right bc
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(negAlternatingL1);
				mesh_[level][imax][j][k]->transport(negAlternatingL1);
			} //!end k
		} //!end j
	} //!end omp

	//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end k
		} //!end j
	} //!end omp

		//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!bottom and top wall
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternatingL1);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	// temporal, spatial interpolation and coupling c2f //
	//--------------------------------------------------//
	//--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
			{
				//!bottom interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

				//!top interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for z = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x y plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				//!interpolation in x y plane
				d2Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for y = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in x z plane
				d2Interpol = 2;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for x = const. //
	//--------------------------------------//

	//--------------------------//
	// mid-nodes in y direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------//
	// mid-nodes in z direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//---------------------------//
	// center-nodes in y z plane //
	//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in y z plane
				d2Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------------------//
	//--------------------------------------------------//
	//        collision and transport (level 1)			//
	//--------------------------------------------------//
	//--------------------------------------------------//

//!collision and transport on level 1 (relevant nodes only)
//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!bottom refined region
				for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
			{
				//!top refined region
				for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!left refined region
				for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
			{
				//!right refined region
				for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end i
	} //!end omp

	//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!left refined region
				for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
			{
				//!right refined region
				for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
				{
					mesh_[level][i][j][k]->collideMRT_RMbasis(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
					mesh_[level][i][j][k]->transport(alternatingL1);
				} //!end i
			} //!end k
		} //!end j
	} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				mesh_[level][i][j][0]->transport(alternatingL1); //!bottom bc
				mesh_[level][i][j][kmax]->transport(alternatingL1); //!top bc
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = 0; i <= imax; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
			{
				mesh_[level][i][0][k]->transport(alternatingL1); //!left bc
				mesh_[level][i][jmax][k]->transport(alternatingL1); //!right bc
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->transport(alternatingL1);
				mesh_[level][imax][j][k]->transport(alternatingL1);
			} //!end k
		} //!end j
	} //!end omp

	//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel;
		//!pressure outlet
		for (int j = 0; j <= jmax; j++) //!all nodes in x direction
		{
			for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
			{
				mesh_[level][0][j][k]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end k
		} //!end j
	} //!end omp

		//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int j = 0; j <= jmax; j++) //!all nodes in y direction
			{
				//!bottom and top wall
				mesh_[level][i][j][0]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do outer loop in parallel
		for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
		{
			for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
			{
				//!left and right wall
				mesh_[level][i][0][k]->LODIpressureOutlet(bc_, alternatingL1);
				mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, alternatingL1);
			} //!end k
		} //!end i
	} //!end omp


	//--------------------------------------------------//
	//--------------------------------------------------//
	//             fine-to-coarse coupling			    //
	//--------------------------------------------------//
	//--------------------------------------------------//

	//!level 0
	level = 0;
	refLayer = bc_.getRefinementLayer();
	refLayerX = bc_.getRefinementLayerX();
	imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
	jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
	kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
			{
				//!bottom interface
				int k = refLayer - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}

				//!top interface
				k = kmax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
		{
			for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}

				//!right interface
				j = jmax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!y direction
		{
			for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!z direction
			{
				//!inlet interface
				int i = refLayerX - overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}

				//!outlet interface
				i = imax - refLayer + overlap;

				if (bc_.getF2Cfilter() == "NOFI")
				{
					mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "LAG")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
				else if (bc_.getF2Cfilter() == "TOU")
				{
					mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
				}
			} //!end k
		} //!end j
	} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             coarse-to-fine coupling			     //
		//--------------------------------------------------//
		//--------------------------------------------------//

	//!level 1
	level = 1;
	refLayer = bc_.getRefinementLayer() * pow(2, level);
	refLayerX = bc_.getRefinementLayerX() * pow(2, level);
	imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
	jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
	kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
			{
				//!bottom interface
				int k = refLayer;
				mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

				//!top interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
			} //!end j
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
			} //!end k
		} //!end i
	} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
			{
				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for z = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
		{
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x y plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				//!interpolation in x y plane
				d2Interpol = 1;

				//!lower interface
				int k = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!upper interface
				k = kmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end j
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for y = const. //
	//--------------------------------------//

		//--------------------------//
		// mid-nodes in x direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
		{
			for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
			{
				//!interpolation in x direction
				d1Interpol = 1;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

		//---------------------------//
		// center-nodes in x z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in x z plane
				d2Interpol = 2;

				//!left interface
				int j = refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!right interface
				j = jmax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
			} //!end k
		} //!end i
	} //!end omp

	//--------------------------------------//
	// spatial interpolation for x = const. //
	//--------------------------------------//

	//--------------------------//
	// mid-nodes in y direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
		{
			for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
			{
				//!interpolation in y direction
				d1Interpol = 2;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//--------------------------//
	// mid-nodes in z direction //
	//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
			{
				//!interpolation in z direction
				d1Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp

	//---------------------------//
	// center-nodes in y z plane //
	//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
	{
#pragma omp for //!do most outer loop in parallel
		for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
		{
			for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
			{
				//!interpolation in y z plane
				d2Interpol = 3;

				//!inlet interface
				int i = refLayerX;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				//!outlet interface
				i = imax - refLayer;
				mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

			} //!end k
		} //!end j
	} //!end omp
} //!end nestedTimeSteppingMRT_RMbasis


void Lattice::nestedTimeSteppingRR(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	if (bc_.getCubicMachCorrection() == "no")
	{
		int d1Interpol; //!switch for 1-D interpolation in different cartesian directions
		int d2Interpol; //!switch for 2-D interpolation in different cartesian directions
		int refLayer; //!number of refined cells
		int refLayerX; //!number of refined cells
		int overlap; //!number of coarse overlap-cells

		//!boolean switches to change distribution arrays
		bool negAlternatingL1 = !alternatingL1;
		bool negAlternatingL0 = !alternatingL0;

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             collision and transport				//
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!current distribution (t=0) level 0: alternating
		//!current distribution (t=0) level 1: alternating

		//!level 0
		int level = 0;
		int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
		refLayer = bc_.getRefinementLayer(); //!number of refined cells
		refLayerX = bc_.getRefinementLayerX(); //!number of refined cells
		overlap = bc_.getOverlap(); //!number of coarse overlap-cells

		/*SAFE PRE-COLLISION DISTRIBUTIONS ON COARSE NODES FOR TEMPORAL INTERPOLATION*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL0, dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!collision and transport on level 0
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!all nodes in z direction
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL0); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL0);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

		//!collision and transport on level 1 (relevant nodes only)
		//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(negAlternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(negAlternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(negAlternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(negAlternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(negAlternatingL1);
					mesh_[level][imax][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

			//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		// temporal, spatial interpolation and coupling c2f //
		//--------------------------------------------------//
		//--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!top interface
					k = kmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!right interface
					j = jmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!outlet interface
					i = imax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in y direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//        collision and transport (level 1)			//
		//--------------------------------------------------//
		//--------------------------------------------------//

//!collision and transport on level 1 (relevant nodes only)
//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(alternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(alternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(alternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(alternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(alternatingL1);
					mesh_[level][imax][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

			//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             fine-to-coarse coupling			    //
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!level 0
		level = 0;
		refLayer = bc_.getRefinementLayer();
		refLayerX = bc_.getRefinementLayerX();
		imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
				{
					//!bottom interface
					int k = refLayer - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!top interface
					k = kmax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!right interface
					j = jmax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!y direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!z direction
				{
					//!inlet interface
					int i = refLayerX - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!outlet interface
					i = imax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end k
			} //!end j
		} //!end omp

			//--------------------------------------------------//
			//--------------------------------------------------//
			//             coarse-to-fine coupling			     //
			//--------------------------------------------------//
			//--------------------------------------------------//

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!top interface
					k = kmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!right interface
					j = jmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!outlet interface
					i = imax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in y direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp
	}
	/* RECURSIVE-REGULARIZED COLLISION MODEL WITH CUBIC MACH CORRECTION */
	else
	{
		int d1Interpol; //!switch for 1-D interpolation in different cartesian directions
		int d2Interpol; //!switch for 2-D interpolation in different cartesian directions
		int refLayer; //!number of refined cells
		int refLayerX; //!number of refined cells
		int overlap; //!number of coarse overlap-cells

		//!boolean switches to change distribution arrays
		bool negAlternatingL1 = !alternatingL1;
		bool negAlternatingL0 = !alternatingL0;

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             collision and transport				//
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!current distribution (t=0) level 0: alternating
		//!current distribution (t=0) level 1: alternating

		//!level 0
		int level = 0;
		int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
		refLayer = bc_.getRefinementLayer(); //!number of refined cells
		refLayerX = bc_.getRefinementLayerX(); //!number of refined cells
		overlap = bc_.getOverlap(); //!number of coarse overlap-cells

		/*SAFE PRE-COLLISION DISTRIBUTIONS ON COARSE NODES FOR TEMPORAL INTERPOLATION*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!z direction
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL0, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!collision and transport on level 0
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!z direction
					{
						//!collide and stream
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL0); //!perform collide at node with hybrid recursive-regularized collision operator (HRR)
						mesh_[level][i][j][k]->transport(negAlternatingL0);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer; k++) //!all fine nodes in bottom interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax; k++) //!all fine nodes in top interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer; j++) //!all fine nodes in left interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax; j++) //!all fine nodes in right interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX; i++) //!all fine nodes at front interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer; i <= imax; i++) //!all fine nodes at rear interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

//!collision and transport on level 1 (relevant nodes only)
//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(negAlternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(negAlternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(negAlternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(negAlternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(negAlternatingL1);
					mesh_[level][imax][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

			//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		// temporal, spatial interpolation and coupling c2f //
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!calculate strain-rate-tensor on coarse partners @ t + 1/2 for temporal interpolation
		level = 0;
		refLayer = bc_.getRefinementLayer();
		refLayerX = bc_.getRefinementLayerX();
		imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					int k = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);

					k = kmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int j = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);

					j = jmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int i = refLayerX;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);

					i = imax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);
				} //!end k
			} //!end j
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!top interface
					k = kmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!right interface
					j = jmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!outlet interface
					i = imax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end j
		} //!end omp


		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in y direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//        collision and transport (level 1)			//
		//--------------------------------------------------//
		//--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer; k++) //!all fine nodes in bottom interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax; k++) //!all fine nodes in top interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer; j++) //!all fine nodes in left interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax; j++) //!all fine nodes in right interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX; i++) //!all fine nodes in inlet interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer; i <= imax; i++) //!all fine nodes in outlet interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!collision and transport on level 1 (relevant nodes only)
		//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
					{
						mesh_[level][i][j][k]->collideRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(alternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(alternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(alternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(alternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(alternatingL1);
					mesh_[level][imax][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

			//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             fine-to-coarse coupling			    //
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!level 0
		level = 0;
		refLayer = bc_.getRefinementLayer();
		refLayerX = bc_.getRefinementLayerX();
		imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

		//!update strain-rate-tensor on fine partners before fine-to-coarse coupling procedure
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
				{
					//!bottom interface
					int k = refLayer - overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }

					//!top interface
					k = kmax - refLayer + overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer - overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }

					//!right interface
					j = jmax - refLayer + overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX - overlap;
					mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);

					//!outlet interface
					i = imax - refLayer + overlap;
					mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

		//!fine-to-coarse coupling
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
				{
					//!bottom interface
					int k = refLayer - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!top interface
					k = kmax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!right interface
					j = jmax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!y direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!z direction
				{
					//!inlet interface
					int i = refLayerX - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!outlet interface
					i = imax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end k
			} //!end j
		} //!end omp


			//--------------------------------------------------//
			//--------------------------------------------------//
			//             coarse-to-fine coupling			     //
			//--------------------------------------------------//
			//--------------------------------------------------//

		//!calculate strain-rate-tensor on coarse partners @ t + 1 for c2f coupling
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					int k = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);

					k = kmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int j = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);

					j = jmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int i = refLayerX;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);

					i = imax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);
				} //!end k
			} //!end j
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!top interface
					k = kmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!right interface
					j = jmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!outlet interface
					i = imax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in y direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp
	} //!end Mach correction
} //!end nestedTimeSteppingRR


//!OUR ALTERNATIVE ALGORITHM WITHOUT FICTITIOUS STREAMING
void Lattice::nestedTimeSteppingHRR(bool& alternatingL1, bool& alternatingL0, Control& bc_) //!perform nested time step with sequence of collide stream step for all level
{
	if (bc_.getCouplingProcedure() == "STD")
	{
		int d1Interpol; //!switch for 1-D interpolation in different cartesian directions
		int d2Interpol; //!switch for 2-D interpolation in different cartesian directions
		int refLayer; //!number of refined cells
		int refLayerX; //!number of refined cells
		int overlap; //!number of coarse overlap-cells

		//!boolean switches to change distribution arrays
		bool negAlternatingL1 = !alternatingL1;
		bool negAlternatingL0 = !alternatingL0;

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             collision and transport				//
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!current distribution (t=0) level 0: alternating
		//!current distribution (t=0) level 1: alternating

		//!level 0
		int level = 0;
		int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
		refLayer = bc_.getRefinementLayer(); //!number of refined cells
		refLayerX = bc_.getRefinementLayerX(); //!number of refined cells
		overlap = bc_.getOverlap(); //!number of coarse overlap-cells

		/*SAFE PRE-COLLISION DISTRIBUTIONS ON COARSE NODES FOR TEMPORAL INTERPOLATION*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!z direction
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL0, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!collision and transport on level 0
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!z direction
					{
						//!collide and stream
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL0); //!perform collide at node with hybrid recursive-regularized collision operator (HRR)
						mesh_[level][i][j][k]->transport(negAlternatingL0);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

		int jmid = jmax / 2;
		int kmid = kmax / 2;

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer; k++) //!all fine nodes in bottom interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax; k++) //!all fine nodes in top interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer; j++) //!all fine nodes in left interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax; j++) //!all fine nodes in right interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX; i++) //!all fine nodes at front interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer; i <= imax; i++) //!all fine nodes at rear interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

	//!collision and transport on level 1 (relevant nodes only)
	//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(negAlternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(negAlternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(negAlternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(negAlternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(negAlternatingL1);
					mesh_[level][imax][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

			//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		// temporal, spatial interpolation and coupling c2f //
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!calculate strain-rate-tensor on coarse partners @ t + 1/2 for temporal interpolation
		level = 0;
		refLayer = bc_.getRefinementLayer();
		refLayerX = bc_.getRefinementLayerX();
		imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					int k = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);

					k = kmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int j = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);

					j = jmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int i = refLayerX;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);

					i = imax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(negAlternatingL0, negAlternatingL1, bc_);
				} //!end k
			} //!end j
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!top interface
					k = kmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!right interface
					j = jmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!outlet interface
					i = imax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFine(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->interpolTimeScalingCoarseToFineMRT_GSbasis(negAlternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in y direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//        collision and transport (level 1)			//
		//--------------------------------------------------//
		//--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer; k++) //!all fine nodes in bottom interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax; k++) //!all fine nodes in top interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer; j++) //!all fine nodes in left interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax; j++) //!all fine nodes in right interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX; i++) //!all fine nodes in inlet interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer; i <= imax; i++) //!all fine nodes in outlet interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

	//!collision and transport on level 1 (relevant nodes only)
	//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(alternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(alternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(alternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(alternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(alternatingL1);
					mesh_[level][imax][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

			//!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             fine-to-coarse coupling			    //
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!level 0
		level = 0;
		refLayer = bc_.getRefinementLayer();
		refLayerX = bc_.getRefinementLayerX();
		imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level

		//!update strain-rate-tensor on fine partners before fine-to-coarse coupling procedure
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
				{
					//!bottom interface
					int k = refLayer - overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }

					//!top interface
					k = kmax - refLayer + overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer - overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }

					//!right interface
					j = jmax - refLayer + overlap;
					if (mesh_[level][i][j][k]->getPartner() != NULL)
					{
						mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
					}
					else { cerr << "No partner!"; }
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX - overlap;
					mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);

					//!outlet interface
					i = imax - refLayer + overlap;
					mesh_[level][i][j][k]->getPartner()->calcStrainRateTensorFDandCubicMachCorrectionFinePartners(alternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

		//!fine-to-coarse coupling
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!all nodes in y direction
				{
					//!bottom interface
					int k = refLayer - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!top interface
					k = kmax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!all nodes in x direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!all nodes in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!right interface
					j = jmax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer - overlap + 1; j <= jmax - refLayer + overlap - 1; j++) //!y direction
			{
				for (int k = refLayer - overlap + 1; k <= kmax - refLayer + overlap - 1; k++) //!z direction
				{
					//!inlet interface
					int i = refLayerX - overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}

					//!outlet interface
					i = imax - refLayer + overlap;

					if (bc_.getF2Cfilter() == "NOFI")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarse(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "LAG")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseLAG(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_LAG_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
					else if (bc_.getF2Cfilter() == "TOU")
					{
						if (bc_.getScalingStyle() == "SRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseTOU(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
						else if (bc_.getScalingStyle() == "MRTSTYLE")
						{
							mesh_[level][i][j][k]->scalingFineToCoarseMRT_TOU_GSbasis(negAlternatingL0, alternatingL1, i, j, k, level, mesh_, bc_);
						}
					}
				} //!end k
			} //!end j
		} //!end omp

			//--------------------------------------------------//
			//--------------------------------------------------//
			//             coarse-to-fine coupling			     //
			//--------------------------------------------------//
			//--------------------------------------------------//

		//!calculate strain-rate-tensor on coarse partners @ t + 1 for c2f coupling
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					int k = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);

					k = kmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int j = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);

					j = jmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int i = refLayerX;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);

					i = imax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(negAlternatingL0, alternatingL1, bc_);
				} //!end k
			} //!end j
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!top interface
					k = kmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!right interface
					j = jmax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

					//!outlet interface
					i = imax - refLayer;
					if (bc_.getScalingStyle() == "SRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFine(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}
					else if (bc_.getScalingStyle() == "MRTSTYLE")
					{
						mesh_[level][i][j][k]->scalingCoarseToFineMRT_GSbasis(alternatingL1, negAlternatingL0, i, j, k, level, mesh_, bc_);
					}

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in y direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp
	} //!end STD coupling


	/*DC COUPLING*/
	else
	{
		int d1Interpol; //!switch for 1-D interpolation in different cartesian directions
		int d2Interpol; //!switch for 2-D interpolation in different cartesian directions
		int refLayer; //!number of refined cells
		int refLayerX; //!number of refined cells
		int overlap; //!number of coarse overlap-cells

		//!boolean switches to change distribution arrays
		bool negAlternatingL1 = !alternatingL1;
		bool negAlternatingL0 = !alternatingL0;

		//--------------------------------------------------//
		//--------------------------------------------------//
		//             collision and transport				//
		//--------------------------------------------------//
		//--------------------------------------------------//

		//!current distribution (t=0) level 0: alternating
		//!current distribution (t=0) level 1: alternating

		//!level 0
		int level = 0;
		int imax = bc_.getNumberOfNodesX() - 1; //!index of last node in x direction in ->this level
		int jmax = bc_.getNumberOfNodesY() - 1; //!index of last node in y direction in ->this level
		int kmax = bc_.getNumberOfNodesZ() - 1; //!index of last node in z direction in ->this level
		refLayer = bc_.getRefinementLayer(); //!number of refined cells
		refLayerX = bc_.getRefinementLayerX(); //!number of refined cells
		overlap = bc_.getOverlap(); //!number of coarse overlap-cells

		/*SAFE PRE-COLLISION DISTRIBUTIONS ON COARSE NODES FOR TEMPORAL INTERPOLATION*/
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!z direction
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL0, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!collision and transport on level 0
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX - overlap; i <= imax - refLayer + overlap; i++) //!x direction
			{
				for (int j = refLayer - overlap; j <= jmax - refLayer + overlap; j++) //!y direction
				{
					for (int k = refLayer - overlap; k <= kmax - refLayer + overlap; k++) //!z direction
					{
						//!collide and stream
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL0); //!perform collide at node with hybrid recursive-regularized collision operator (HRR)
						mesh_[level][i][j][k]->transport(negAlternatingL0);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer; k++) //!all fine nodes in bottom interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax; k++) //!all fine nodes in top interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer; j++) //!all fine nodes in left interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax; j++) //!all fine nodes in right interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX; i++) //!all fine nodes at front interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer; i <= imax; i++) //!all fine nodes at rear interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(alternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

	//!collision and transport on level 1 (relevant nodes only)
	//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at outlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, alternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(negAlternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

        //!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(negAlternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(negAlternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(negAlternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(negAlternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(negAlternatingL1);
					mesh_[level][imax][j][k]->transport(negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		 //!perform bounce-back on level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, negAlternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, negAlternatingL1);
				} //!end k
			} //!end i
		} //!end omp

		//-----------------------------------------------------//
		//-----------------------------------------------------//
		// temporal, spatial interpolation and direct coupling //
		//-----------------------------------------------------//
		//-----------------------------------------------------//

		#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolTimeDirectCoupling(negAlternatingL1, negAlternatingL0, bc_);

					//!top interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolTimeDirectCoupling(negAlternatingL1, negAlternatingL0, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolTimeDirectCoupling(negAlternatingL1, negAlternatingL0, bc_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolTimeDirectCoupling(negAlternatingL1, negAlternatingL0, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolTimeDirectCoupling(negAlternatingL1, negAlternatingL0, bc_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolTimeDirectCoupling(negAlternatingL1, negAlternatingL0, bc_);

				} //!end k
			} //!end j
		} //!end omp

		//!calculate strain-rate-tensor on fine interface nodes @ t + 1/2
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all nodes in y direction
				{
					int k = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, negAlternatingL1, bc_);

					k = kmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, negAlternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all nodes in z direction
				{
					int j = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, negAlternatingL1, bc_);

					j = jmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, negAlternatingL1, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all nodes in y direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all nodes in z direction
				{
					int i = refLayerX;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, negAlternatingL1, bc_);

					i = imax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, negAlternatingL1, bc_);
				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

            //--------------------------//
            // mid-nodes in y direction //
            //--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(negAlternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(negAlternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------------------------------//
		//--------------------------------------------------//
		//        collision and transport (level 1)			//
		//--------------------------------------------------//
		//--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!bottom refined region
					for (int k = 0; k <= refLayer; k++) //!all fine nodes in bottom interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				for (int j = 0; j <= jmax; j++) //!y direction
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax; k++) //!all fine nodes in top interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
						}
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int j = 0; j <= refLayer; j++) //!all fine nodes in left interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!x direction
			{
				//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax; j++) //!all fine nodes in right interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end j
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!inlet refined region
					for (int i = 0; i <= refLayerX; i++) //!all fine nodes in inlet interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			//!inlet and outlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!outlet refined region
					for (int i = imax - refLayer; i <= imax; i++) //!all fine nodes in outlet interface
					{
						for (int dir = 0; dir < 19; dir++)
						{
							double distrTemp = mesh_[level][i][j][k]->getDistribution(negAlternatingL1, dir); //!set temporary distribution
							mesh_[level][i][j][k]->setDistributionPreCol(dir, distrTemp);
							double psiTemp = mesh_[level][i][j][k]->getCubicMachCorrection(dir); //!initialize temporary distribution
							mesh_[level][i][j][k]->setCubicMachCorrectionPreCol(dir, psiTemp);
						}
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

	//!collision and transport on level 1 (relevant nodes only)
	//!bottom and top refined region
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!bottom refined region
					for (int k = 1; k <= refLayer; k++) //!all fine nodes at bottom interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 1; j <= jmax - 1; j++) //!all fine nodes in y direction except wall nodes
				{
					//!top refined region
					for (int k = kmax - refLayer; k <= kmax - 1; k++) //!all fine nodes at top interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!left and right refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!left refined region
					for (int j = 1; j <= refLayer; j++) //!all fine nodes at left interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 1; i <= imax -1; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all fine nodes in z direction except wall nodes and nodes already considered before
				{
					//!right refined region
					for (int j = jmax - refLayer; j <= jmax - 1; j++) //!all fine nodes at right interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end k
				} //!end j
			} //!end i
		} //!end omp

		//!inlet refined region (KEEP IN MIND: SOME NODES ALREADY CONSIDERED IN INTERFACE REGIONS ABOVE!)
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = 1; i <= refLayerX; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!z direction
				{
					//!left refined region
					for (int i = imax - refLayer; i <= imax - 1; i++) //!all fine nodes at inlet interface
					{
						mesh_[level][i][j][k]->collideHRR(bc_, negAlternatingL1); //!perform collide at node with BGK collision operator
						mesh_[level][i][j][k]->transport(alternatingL1);
					} //!end i
				} //!end k
			} //!end j
		} //!end omp

		//!let boundary nodes first perform transport step to update missing populations in neighbors
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					mesh_[level][i][j][0]->transport(alternatingL1); //!bottom bc
					mesh_[level][i][j][kmax]->transport(alternatingL1); //!top bc
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = 0; i <= imax; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction except boundary
				{
					mesh_[level][i][0][k]->transport(alternatingL1); //!left bc
					mesh_[level][i][jmax][k]->transport(alternatingL1); //!right bc
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			for (int j = 1; j <= jmax - 1; j++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->transport(alternatingL1);
					mesh_[level][imax][j][k]->transport(alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

		//!pressure outlet in level 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel;
			//!pressure outlet
			for (int j = 0; j <= jmax; j++) //!all nodes in x direction
			{
				for (int k = 0; k <= kmax; k++)  //!all nodes in y direction
				{
					mesh_[level][0][j][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][imax][j][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end j
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int j = 0; j <= jmax; j++) //!all nodes in y direction
				{
					//!bottom and top wall
					mesh_[level][i][j][0]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][j][kmax]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = 1; i <= imax - 1; i++) //!all nodes in x direction
			{
				for (int k = 1; k <= kmax - 1; k++) //!all nodes in z direction
				{
					//!left and right wall
					mesh_[level][i][0][k]->LODIpressureOutlet(bc_, alternatingL1);
					mesh_[level][i][jmax][k]->LODIpressureOutlet(bc_, alternatingL1);
				} //!end k
			} //!end i
		} //!end omp

        //--------------------------------------------------//
        //--------------------------------------------------//
        //					direct coupling			        //
        //--------------------------------------------------//
        //--------------------------------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!all fine nodes at coarse positions in y direction
				{
					//!bottom interface
					int k = refLayer;
					mesh_[level][i][j][k]->directCoupling(alternatingL1, negAlternatingL0, bc_);

					//!top interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->directCoupling(alternatingL1, negAlternatingL0, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!all fine nodes at coarse positions in x direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->directCoupling(alternatingL1, negAlternatingL0, bc_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->directCoupling(alternatingL1, negAlternatingL0, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 2; j <= jmax - refLayer - 2; j += 2) //!all fine nodes at coarse positions in y direction except inlet and outlet
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!all fine nodes at coarse positions in z direction except nodes already considered before
				{
					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->directCoupling(alternatingL1, negAlternatingL0, bc_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->directCoupling(alternatingL1, negAlternatingL0, bc_);

				} //!end k
			} //!end j
		} //!end omp

		//!level 0
		level = 0;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

		//!calculate strain-rate-tensor on coarse partners @ t + 1
#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j++) //!all nodes in y direction
				{
					int k = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, alternatingL1, bc_);

					k = kmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, alternatingL1, bc_);
				} //!end j
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i++) //!all nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int j = refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, alternatingL1, bc_);

					j = jmax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, alternatingL1, bc_);
				} //!end k
			} //!end i
		} //!end omp

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j++) //!all nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k++) //!all nodes in z direction
				{
					int i = refLayerX;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, alternatingL1, bc_);

					i = imax - refLayer;
					mesh_[level][i][j][k]->calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(negAlternatingL0, alternatingL1, bc_);
				} //!end k
			} //!end j
		} //!end omp

		//!level 1
		level = 1;
		refLayer = bc_.getRefinementLayer() * pow(2, level);
		refLayerX = bc_.getRefinementLayerX() * pow(2, level);
		imax = (bc_.getNumberOfNodesX() - 1) * pow(2., level); //!index of last node in x direction in ->this level
		jmax = (bc_.getNumberOfNodesY() - 1) * pow(2., level); //!index of last node in y direction in ->this level
		kmax = (bc_.getNumberOfNodesZ() - 1) * pow(2., level); //!index of last node in z direction in ->this level

		//--------------------------------------//
		// spatial interpolation for z = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all bottom interface mid-nodes in x direction
			{
				for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in y direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all mid-nodes in y direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x y plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
				{
					//!interpolation in x y plane
					d2Interpol = 1;

					//!lower interface
					int k = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!upper interface
					k = kmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end j
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for y = const. //
		//--------------------------------------//

			//--------------------------//
			// mid-nodes in x direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all left interface mid-nodes in x direction
			{
				for (int k = refLayer + 2; k <= kmax - refLayer - 2; k += 2) //!z direction
				{
					//!interpolation in x direction
					d1Interpol = 1;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//--------------------------//
			// mid-nodes in z direction //
			//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX; i <= imax - refLayer; i += 2) //!x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

			//---------------------------//
			// center-nodes in x z plane //
			//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int i = refLayerX + 1; i <= imax - refLayer - 1; i += 2) //!all center-nodes in x direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in x z plane
					d2Interpol = 2;

					//!left interface
					int j = refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!right interface
					j = jmax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);
				} //!end k
			} //!end i
		} //!end omp

		//--------------------------------------//
		// spatial interpolation for x = const. //
		//--------------------------------------//

		//--------------------------//
		// mid-nodes in y direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all inlet interface mid-nodes in y direction
			{
				for (int k = refLayer; k <= kmax - refLayer; k += 2) //!z direction
				{
					//!interpolation in y direction
					d1Interpol = 2;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//--------------------------//
		// mid-nodes in z direction //
		//--------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer; j <= jmax - refLayer; j += 2) //!y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all mid-nodes in z direction
				{
					//!interpolation in z direction
					d1Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceMidNode(alternatingL1, i, j, k, d1Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp

		//---------------------------//
		// center-nodes in y z plane //
		//---------------------------//

#pragma omp parallel num_threads(bc_.getNumberOfThreads()) //!shared memory parallelization
		{
#pragma omp for //!do most outer loop in parallel
			for (int j = refLayer + 1; j <= jmax - refLayer - 1; j += 2) //!all center-nodes in y direction
			{
				for (int k = refLayer + 1; k <= kmax - refLayer - 1; k += 2) //!all center-nodes in z direction
				{
					//!interpolation in y z plane
					d2Interpol = 3;

					//!inlet interface
					int i = refLayerX;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

					//!outlet interface
					i = imax - refLayer;
					mesh_[level][i][j][k]->interpolSpaceCenterNode(alternatingL1, i, j, k, d2Interpol, level, bc_, mesh_);

				} //!end k
			} //!end j
		} //!end omp
	} //!end direct coupling
} //!end nestedTimeSteppingHRR


double Lattice::readValue(string filename, string keyword) //!filename: name of initialization file; keyword: string to be read from file; return value: (double) to corresponding keyword
{
	ifstream file;
	string s;
	file.open(filename.c_str());
	istringstream sin;
	string op;
	double number = 0.;

	if (!file)
	{
		cerr << "Error reading initialization file " << filename << "!";
		exit(1);
	}

	while (file.good())
	{
		op = " ";
		getline(file, s);
		sin.clear();
		sin.str(s);
		sin >> op;
		if (op == keyword)
		{
			sin >> number;
			return number;
		}
	}
	cerr << "No Value for keyword \"" + keyword + "\" found in " << filename << "!";
	exit(1);
	return number;
}

