#include "Node.h"
#include "Vector3D.H"
#include <vector>
#include <map>
#include <cmath>

///////////////////////////////////////////////////////////////////
//!MACROSCOPIC VALUES
///////////////////////////////////////////////////////////////////

double Node::calcDensity(Control& bc, bool& alternating)
{
	double dens = 0;

	if (bc.getCubicMachCorrection() == "no")
	{
		for (int dir = 0; dir < 19; dir++)
		{
			dens += distr_[alternating][dir];
		}
	}
	else
	{
		for (int dir = 0; dir < 19; dir++)
		{
			dens += distr_[alternating][dir] + 0.5 * timeStep_ * psi_[dir];
		}
	}
	return dens;
}

double Node::calcDensityPreCol(Control& bc)
{
	double dens = 0;

	if (bc.getCubicMachCorrection() == "no")
	{
		for (int dir = 0; dir < 19; dir++)
		{
			dens += distrPreCol_[dir];
		}
	}
	else
	{
		for (int dir = 0; dir < 19; dir++)
		{
			dens += distrPreCol_[dir] + 0.5 * timeStep_ * psiPreCol_[dir];
		}
	}
	return dens;
}

double Node::calcDensityFictitious()
{
	double dens = 0;

	for (int dir = 0; dir < 19; dir++)
	{
		dens += distrFictitious_[dir];
	}
	return dens;
}

double Node::calcPressure(bool& alternating, Control& bc)
{
	double dens = 0;

	if (bc.getCubicMachCorrection() == "no")
	{
		for (int dir = 0; dir < 19; dir++)
		{
			dens += distr_[alternating][dir];
		}
	}
	else
	{
		for (int dir = 0; dir < 19; dir++)
		{
			dens += distr_[alternating][dir] + 0.5 * timeStep_ * psi_[dir];
		}
	}
	return dens * (pow(bc.getSpeedOfSound(), 2.));
}

void Node::calcVelocity(Control& bc, bool& alternating, double& dens, double* velo)   //!choose distribution by alternating
{
	velo[0] = 0.; velo[1] = 0.; velo[2] = 0.; //!init variable

	//!get lattice velocities
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();

	if (bc.getCubicMachCorrection() == "no")
	{
		for (int dir = 0; dir < 19; dir++) //!loop over all directions (D3Q19)
		{
			velo[0] += xsiX[dir] * distr_[alternating][dir];
			velo[1] += xsiY[dir] * distr_[alternating][dir];
			velo[2] += xsiZ[dir] * distr_[alternating][dir];
		}
	}
	else
	{
		for (int dir = 0; dir < 19; dir++) //!loop over all directions (D3Q19)
		{
			velo[0] += xsiX[dir] * distr_[alternating][dir] + 0.5 * timeStep_ * xsiX[dir] * psi_[dir];
			velo[1] += xsiY[dir] * distr_[alternating][dir] + 0.5 * timeStep_ * xsiY[dir] * psi_[dir];
			velo[2] += xsiZ[dir] * distr_[alternating][dir] + 0.5 * timeStep_ * xsiZ[dir] * psi_[dir];
		}
	}

	velo[0] /= dens;
	velo[1] /= dens;
	velo[2] /= dens;

	return; //!result is passed by reference velo
}


void Node::calcVelocityPreCol(Control& bc, double& dens, double* velo)   //!choose distribution by alternating
{
	velo[0] = 0.; velo[1] = 0.; velo[2] = 0.; //!init variable

	//!get lattice velocities
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();

	if (bc.getCubicMachCorrection() == "no")
	{
		for (int dir = 0; dir < 19; dir++) //!loop over all directions (D3Q19)
		{
			velo[0] += xsiX[dir] * distrPreCol_[dir];
			velo[1] += xsiY[dir] * distrPreCol_[dir];
			velo[2] += xsiZ[dir] * distrPreCol_[dir];
		}
	}
	else
	{
		for (int dir = 0; dir < 19; dir++) //!loop over all directions (D3Q19)
		{
			velo[0] += xsiX[dir] * distrPreCol_[dir] + 0.5 * timeStep_ * xsiX[dir] * psiPreCol_[dir];
			velo[1] += xsiY[dir] * distrPreCol_[dir] + 0.5 * timeStep_ * xsiY[dir] * psiPreCol_[dir];
			velo[2] += xsiZ[dir] * distrPreCol_[dir] + 0.5 * timeStep_ * xsiZ[dir] * psiPreCol_[dir];
		}
	}

	velo[0] /= dens;
	velo[1] /= dens;
	velo[2] /= dens;

	return; //!result is passed by reference velo
}

void Node::calcVelocityFictitious(Control& bc, double& dens, double* velo)   //!choose distribution by alternating
{
	velo[0] = 0.; velo[1] = 0.; velo[2] = 0.; //!init variable

	//!get lattice velocities
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();

	for (int dir = 0; dir < 19; dir++) //!loop over all directions (D3Q19)
	{
		velo[0] += xsiX[dir] * distrFictitious_[dir];
		velo[1] += xsiY[dir] * distrFictitious_[dir];
		velo[2] += xsiZ[dir] * distrFictitious_[dir];
	}

	velo[0] /= dens;
	velo[1] /= dens;
	velo[2] /= dens;

	return; //!result is passed by reference velo
}

void Node::pressureOutlet(Control& bc, bool& alternating)
{
	bool negAlternating = !alternating;

	//!read prescribed inlet density
	double dens = bc.getDensity();

	//!interpolate density from next fluid node neighbor in negative x direction
	double velo[3];
	getNeighbor(2)->calcVelocity(bc, negAlternating, dens, velo);

	//!calculate equilibrium distribution from prescribed inlet velocity and interpolated density
	double equil[19];
    if (bc.getOrderOfEquilibrium() == 2)
	{
        this->calcEquilibrium2_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 3)
	{
        this->calcEquilibrium3_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 4)
	{
        this->calcEquilibrium4_GHbasis(equil, bc, dens, velo);
	}
	else
	{
        cerr << "Wrong order of equilibrium!" << flush; exit(1);
	}

	////!calculate deviatoric stress tensor (deviatoric stress tensor is a member of node class)
	//this->calcDevStressTensor(alternating, equil, bc, velo);

	////!calculate first-order off-equilibrium distribution from recursive-regularized approach (third-order hermite polynomial basis)
	//double nEquilRR[19];
	//this->calcRRoffEquilibrium(alternating, nEquilRR, equil, bc, velo);

	for (int dir = 0; dir < 18; dir++)
	{
		this->setDistribution(alternating, dir, equil[dir]);
		this->setDistribution(negAlternating, dir, equil[dir]);
		if (neighbors_[dir] != NULL)
		{
			//neighbors_[dir]->setDistribution(alternating, dir, equil[dir] + nEquilRR[dir]);
			neighbors_[dir]->setDistribution(alternating, dir, equil[dir]);
		}
	}
}


void Node::LODIpressureOutlet(Control& bc, bool& alternating)
{
	//!solve LODI equations to achieve non reflecting BC; see Feuchter: Direct aeroacoustic simulation with a cumulant Lattice-Boltzmann model. Computers and Fluids 224 (2021) 104970, p. 8
	//!here for simplicity only pressure outlet with straight normal directions are considered! (i.e. flow in +/- x,y,z-direction!)
	bool negAlternating = !alternating;
	double density = 0.;
	double velocity[3] = { 0. };
	static double speedOfSound = bc.getSpeedOfSound();
	static double speedOfSoundSq = pow(speedOfSound, 2.);
	static double spacing = bc.getSpacing();
	double localSpacing = spacing * pow(2., -level_);

	//!initialize transverse gradients
	double densityT11 = 0.; double densityT12 = 0.; double densityT21 = 0.; double densityT22 = 0.;
	double densityGradientT[2] = { 0. };
	double velocityT11[3] = { 0. };
	double velocityT12[3] = { 0. };
	double velocityT21[3] = { 0. };
	double velocityT22[3] = { 0. };
	double velocityGradientT[2][3] = { 0. };

	//!calculate amplitude of outgoing acoustic wave
	double densityBoundaryNode = this->calcDensity(bc, negAlternating);
	double veloBoundaryNode[3] = { 0. };
	this->calcVelocity(bc, negAlternating, densityBoundaryNode, veloBoundaryNode); //!use calcVelocity of Node to get current value not boundary condition value itself
	double densityGradient = 0.; //!gradient of density normal to boundary
	double velocityGradient[3] = { 0. }; //!gradient of x-velo, y-velo and z-velo normal to boundary

	if (this->getNormalNeighbor() != NULL)
	{
		//!calculate density gradient
		double densityNeighbor = this->getNormalNeighbor()->calcDensity(bc, alternating);
		double density2ndNeighbor = this->getNormalNeighborSnd()->calcDensity(bc, alternating);
		densityGradient = (4. * densityNeighbor - density2ndNeighbor - 3. * densityBoundaryNode) / (-2. * localSpacing); //!calculate density gradient by second order one sided finite difference (left-sided FD) --> consider - in denominator due to backward difference --> works only for boundaries with pos. normal vector

		//!calculate velocity gradient
		double velocityNeighbor[3] = { 0. };
		this->getNormalNeighbor()->calcVelocity(bc, alternating, densityNeighbor, velocityNeighbor);
		double velocity2ndNeighbor[3] = { 0. };
		this->getNormalNeighborSnd()->calcVelocity(bc, alternating, density2ndNeighbor, velocity2ndNeighbor);
		velocityGradient[0] = (4. * velocityNeighbor[0] - velocity2ndNeighbor[0] - 3. * veloBoundaryNode[0]) / (-2. * localSpacing); //!calculate x-velo gradient by second order one sided finite difference (left-sided FD) --> consider - in denominator due to backward difference --> works only for boundaries with pos. normal vector; for neg. normal vector *-1
		velocityGradient[1] = (4. * velocityNeighbor[1] - velocity2ndNeighbor[1] - 3. * veloBoundaryNode[1]) / (-2. * localSpacing); //!calculate y-velo gradient by second order one sided finite difference (left-sided FD)
		velocityGradient[2] = (4. * velocityNeighbor[2] - velocity2ndNeighbor[2] - 3. * veloBoundaryNode[2]) / (-2. * localSpacing); //!calculate z-velo gradient by second order one sided finite difference (left-sided FD)

		//!calculate transverse gradients (CENTRAL DIFFERENCE SCHEME)
		//!tangent neighbors refer to all four cartesian in-plane neighbors of this->node (corners/edges not yet included!)
		//!for tangent neighbors x comes before y comes before z
		if ((!this->getIsCorner()) && (bc.getTransverseInclusion() == "yes"))
		{
			densityT11 = this->getTangentNeighbor1(0)->calcDensity(bc, negAlternating);
			densityT12 = this->getTangentNeighbor1(1)->calcDensity(bc, negAlternating);
			densityT21 = this->getTangentNeighbor2(0)->calcDensity(bc, negAlternating);
			densityT22 = this->getTangentNeighbor2(1)->calcDensity(bc, negAlternating);

			this->getTangentNeighbor1(0)->calcVelocity(bc, negAlternating, densityT11, velocityT11);
			this->getTangentNeighbor1(1)->calcVelocity(bc, negAlternating, densityT12, velocityT12);
			this->getTangentNeighbor2(0)->calcVelocity(bc, negAlternating, densityT21, velocityT21);
			this->getTangentNeighbor2(1)->calcVelocity(bc, negAlternating, densityT22, velocityT22);

			densityGradientT[0] = (densityT11 - densityT12) / (2. * localSpacing);
			densityGradientT[1] = (densityT21 - densityT22) / (2. * localSpacing);
			velocityGradientT[0][0] = (velocityT11[0] - velocityT12[0]) / (2. * localSpacing);
			velocityGradientT[0][1] = (velocityT11[1] - velocityT12[1]) / (2. * localSpacing);
			velocityGradientT[0][2] = (velocityT11[2] - velocityT12[2]) / (2. * localSpacing);
			velocityGradientT[1][0] = (velocityT21[0] - velocityT22[0]) / (2. * localSpacing);
			velocityGradientT[1][1] = (velocityT21[1] - velocityT22[1]) / (2. * localSpacing);
			velocityGradientT[1][2] = (velocityT21[2] - velocityT22[2]) / (2. * localSpacing);
		}
	}
	else
	{
		cerr << "Normal neighbor is missing!" << endl;
	}

	//!currently only boundaries with +/- x,y,z normals are support by LODI
	//!check for normal direction
	if (this->getNormalDir() == 0) //!boundary has -x normal
	{
		densityGradient *= -1.; //!use right sided FD, see comment above
		velocityGradient[0] *= -1.; //!use right sided FD, see comment above
		velocityGradient[1] *= -1.; //!use right sided FD, see comment above
		velocityGradient[2] *= -1.; //!use right sided FD, see comment above
		double domainLength = bc.getNumberOfNodesX() * spacing;
		//!switch velocity directions if necessary, i.e. use local coordinate system to make use of identical calculus, i.e. function
		double switchedVelocityGradient[3] = { velocityGradient[0], velocityGradient[1], velocityGradient[2] };
		double switchedVelocityBoundaryNode[3] = { veloBoundaryNode[0], veloBoundaryNode[1], veloBoundaryNode[2] };
		LODIpressureCalcVelocity(bc, domainLength, density, densityBoundaryNode, densityGradient, switchedVelocityGradient, switchedVelocityBoundaryNode, velocity, true, densityGradientT, velocityGradientT);
		//!switch velocity directions back again, to be compliant to global coordinate system
		double temp[3] = { velocity[0], velocity[1], velocity[2] };
		velocity[0] = temp[0];
		velocity[1] = temp[1];
		velocity[2] = temp[2];
	}
	else if (this->getNormalDir() == 2) //!boundary has +x normal
	{
		double domainLength = bc.getNumberOfNodesX() * spacing;
		LODIpressureCalcVelocity(bc, domainLength, density, densityBoundaryNode, densityGradient, velocityGradient, veloBoundaryNode, velocity, false, densityGradientT, velocityGradientT);
	}
	else if (this->getNormalDir() == 1) //boundary has -y normal
	{
		densityGradient *= -1.; //!use right sided FD, see comment above
		velocityGradient[0] *= -1.; //!use right sided FD, see comment above
		velocityGradient[1] *= -1.; //!use right sided FD, see comment above
		velocityGradient[2] *= -1.; //!use right sided FD, see comment above
		double domainLength = bc.getNumberOfNodesY() * spacing;
		//!switch x and y coordinate, i.e. use local coordinate system to make use of identical calculus, i.e. function
		double switchedVelocityGradient[3] = { velocityGradient[1], velocityGradient[0], velocityGradient[2] };
		double switchedVelocityGradientT[2][3] = { 0. };
		switchedVelocityGradientT[0][0] = velocityGradientT[0][1];
		switchedVelocityGradientT[0][1] = velocityGradientT[0][0];
		switchedVelocityGradientT[0][2] = velocityGradientT[0][2];
		switchedVelocityGradientT[1][0] = velocityGradientT[1][1];
		switchedVelocityGradientT[1][1] = velocityGradientT[1][0];
		switchedVelocityGradientT[1][2] = velocityGradientT[1][2];
		double switchedVeloBoundaryNode[3] = { veloBoundaryNode[1], veloBoundaryNode[0], veloBoundaryNode[2] };
		LODIpressureCalcVelocity(bc, domainLength, density, densityBoundaryNode, densityGradient, switchedVelocityGradient, switchedVeloBoundaryNode, velocity, true, densityGradientT, switchedVelocityGradientT);
		//switch x and y direction again, to be compliant to global coordinate system
		double temp[3] = { velocity[0], velocity[1], velocity[2] };
		velocity[0] = temp[1];
		velocity[1] = temp[0];
		velocity[2] = temp[2];
	}
	else if (this->getNormalDir() == 3) //!boundary has +y normal
	{
		double domainLength = bc.getNumberOfNodesY() * spacing;
		//!switch x and y coordinate, i.e. use local coordinate system to make use of identical calculus, i.e. function
		double switchedVelocityGradient[3] = { velocityGradient[1], velocityGradient[0], velocityGradient[2] };
		double switchedVelocityGradientT[2][3] = { 0. };
		switchedVelocityGradientT[0][0] = velocityGradientT[0][1];
		switchedVelocityGradientT[0][1] = velocityGradientT[0][0];
		switchedVelocityGradientT[0][2] = velocityGradientT[0][2];
		switchedVelocityGradientT[1][0] = velocityGradientT[1][1];
		switchedVelocityGradientT[1][1] = velocityGradientT[1][0];
		switchedVelocityGradientT[1][2] = velocityGradientT[1][2];
		double switchedVeloBoundaryNode[3] = { veloBoundaryNode[1], veloBoundaryNode[0], veloBoundaryNode[2] };
		LODIpressureCalcVelocity(bc, domainLength, density, densityBoundaryNode, densityGradient, switchedVelocityGradient, switchedVeloBoundaryNode, velocity, false, densityGradientT, switchedVelocityGradientT);
		//!switch x and y direction again, to be compliant to global coordinate system
		double temp[3] = { velocity[0], velocity[1], velocity[2] };
		velocity[0] = temp[1];
		velocity[1] = temp[0];
		velocity[2] = temp[2];
	}
	else if (this->getNormalDir() == 4) //!boundary has -z normal
	{
		densityGradient *= -1.; //!use right sided FD, see comment above
		velocityGradient[0] *= -1.; //!use right sided FD, see comment above
		velocityGradient[1] *= -1.; //!use right sided FD, see comment above
		velocityGradient[2] *= -1.; //!use right sided FD, see comment above
		double domainLength = bc.getNumberOfNodesZ() * spacing;
		//!switch x and z as well as z and y coordinate, i.e. use local coordinate system to make use of identical calculus, i.e. function
		double switchedVelocityGradient[3] = { velocityGradient[2], velocityGradient[0], velocityGradient[1] };
		double switchedVelocityGradientT[2][3] = { 0. };
		switchedVelocityGradientT[0][0] = velocityGradientT[0][2];
		switchedVelocityGradientT[0][1] = velocityGradientT[0][0];
		switchedVelocityGradientT[0][2] = velocityGradientT[0][1];
		switchedVelocityGradientT[1][0] = velocityGradientT[1][2];
		switchedVelocityGradientT[1][1] = velocityGradientT[1][0];
		switchedVelocityGradientT[1][2] = velocityGradientT[1][1];
		double switchedVeloBoundaryNode[3] = { veloBoundaryNode[2], veloBoundaryNode[0], veloBoundaryNode[1] };
		LODIpressureCalcVelocity(bc, domainLength, density, densityBoundaryNode, densityGradient, switchedVelocityGradient, switchedVeloBoundaryNode, velocity, true, densityGradientT, switchedVelocityGradientT);
		//!switch components back again, to be compliant to global coordinate system
		double temp[3] = { velocity[0], velocity[1], velocity[2] };
		velocity[0] = temp[1];
		velocity[1] = temp[2];
		velocity[2] = temp[0];
	}
	else if (this->getNormalDir() == 5) //!boundary has +z normal
	{
		double domainLength = bc.getNumberOfNodesZ() * spacing;
		//!switch x and z as well as z and y coordinate, i.e. use local coordinate system to make use of identical calculus, i.e. function
		double switchedVelocityGradient[3] = { velocityGradient[2], velocityGradient[0], velocityGradient[1] };
		double switchedVelocityGradientT[2][3] = { 0. };
		switchedVelocityGradientT[0][0] = velocityGradientT[0][2];
		switchedVelocityGradientT[0][1] = velocityGradientT[0][0];
		switchedVelocityGradientT[0][2] = velocityGradientT[0][1];
		switchedVelocityGradientT[1][0] = velocityGradientT[1][2];
		switchedVelocityGradientT[1][1] = velocityGradientT[1][0];
		switchedVelocityGradientT[1][2] = velocityGradientT[1][1];
		double switchedVeloBoundaryNode[3] = { veloBoundaryNode[2], veloBoundaryNode[0], veloBoundaryNode[1] };
		LODIpressureCalcVelocity(bc, domainLength, density, densityBoundaryNode, densityGradient, switchedVelocityGradient, switchedVeloBoundaryNode, velocity, false, densityGradientT, switchedVelocityGradientT);
		//!switch components back again, to be compliant to global coordinate system
		double temp[3] = { velocity[0], velocity[1], velocity[2] };
		velocity[0] = temp[1];
		velocity[1] = temp[2];
		velocity[2] = temp[0];
	}
	else
	{
		cerr << "Some node has no normalDir!" << endl;
	}

	//!calculate equilibrium distribution from prescribed inlet velocity and interpolated density
	double equil[19];
    if (bc.getOrderOfEquilibrium() == 2)
	{
        this->calcEquilibrium2_GHbasis(equil, bc, density, velocity);
	}
	else if (bc.getOrderOfEquilibrium() == 3)
	{
        this->calcEquilibrium3_GHbasis(equil, bc, density, velocity);
	}
	else if (bc.getOrderOfEquilibrium() == 4)
	{
        this->calcEquilibrium4_GHbasis(equil, bc, density, velocity);
	}
	else
	{
        cerr << "Wrong order of equilibrium!" << flush; exit(1);
	}

	for (int dir = 0; dir < 19; dir++)
	{
		this->setDistribution(alternating, dir, equil[dir]);
	}
}//!end LODIpressureOutlet



void Node::LODIpressureCalcVelocity(Control& bc, double& domainLength, double& density, double& densityBoundaryNode, double& densityGradient, double* velocityGradient, double* veloBoundaryNode, double* velocity, bool negNormal, double* densityGradientT, double velocityGradientT[][3])
{
	static double speedOfSound = bc.getSpeedOfSound();
	static double speedOfSoundSq = pow(speedOfSound, 2.);

	//!calculate amplitudes
	double csDivL = speedOfSound / domainLength;
	static double Ma = bc.getMachNumber(); //!approximated max. Mach number in domain
	static double MaSq = Ma * Ma;

	//!calculate amplitude of incoming wave by relaxation
	static double K1 = bc.getLODIK1relaxation(); //!relaxation rates for incoming longitudal waves
	static double K2 = Ma; //!relaxation rates for incoming transverse waves (should be equal to average MACH-number)
	double densityOutlet = bc.getDensity(); //!calculate density at outlet
	double veloX = 0.;
	double veloY = 0.;
	double veloZ = 0.;
	double timeStep = bc.getTimeStep() / pow(2., level_);

	double L1, L5, T1, T5;
	double Texact = 0.;

	//!calculate transverse Terms T1 and T5
	T1 = -(speedOfSoundSq * (densityBoundaryNode * (velocityGradientT[0][1] + velocityGradientT[1][2]) + veloBoundaryNode[1] * densityGradientT[0] + veloBoundaryNode[2] * densityGradientT[1]) - densityBoundaryNode * speedOfSound * (veloBoundaryNode[1] * velocityGradientT[0][0] + veloBoundaryNode[2] * velocityGradientT[1][0]));
	T5 = -(speedOfSoundSq * (densityBoundaryNode * (velocityGradientT[0][1] + velocityGradientT[1][2]) + veloBoundaryNode[1] * densityGradientT[0] + veloBoundaryNode[2] * densityGradientT[1]) + densityBoundaryNode * speedOfSound * (veloBoundaryNode[1] * velocityGradientT[0][0] + veloBoundaryNode[2] * velocityGradientT[1][0]));

	if (negNormal)
	{
		//!for negative boundary normal L1 and L5 are switched --> L5 is incoming, L1 outgoing
		L5 = K1 * (1 - MaSq) * csDivL * speedOfSoundSq * (densityBoundaryNode - densityOutlet) - K2 * (T5 - Texact) + T5;
		L1 = (veloBoundaryNode[0] - speedOfSound) * speedOfSound * (speedOfSound * densityGradient - densityBoundaryNode * velocityGradient[0]);
	}
	else
	{
		L1 = K1 * (1 - MaSq) * csDivL * speedOfSoundSq * (densityBoundaryNode - densityOutlet) - K2 * (T1 - Texact) + T1;
		//!calculate amplitudes of outgoing waves
		L5 = (veloBoundaryNode[0] + speedOfSound) * speedOfSound * (speedOfSound * densityGradient + densityBoundaryNode * velocityGradient[0]);
	}

	//!calculate amplitudes of remaining waves
	double L3 = veloBoundaryNode[0] * velocityGradient[1];
	double L4 = veloBoundaryNode[0] * velocityGradient[2];

	//!calculate transverse terms T3 and T4
	double T3 = -(veloBoundaryNode[1] * velocityGradientT[0][1] + veloBoundaryNode[2] * velocityGradientT[1][1]) - (speedOfSoundSq / densityBoundaryNode) * densityGradientT[0];
	double T4 = -(veloBoundaryNode[1] * velocityGradientT[0][2] + veloBoundaryNode[2] * velocityGradientT[1][2]) - (speedOfSoundSq / densityBoundaryNode) * densityGradientT[1];

	//!calculate velocity components and density from LODI equations --> use first order finite differences for time derivatives
	density = densityBoundaryNode - timeStep / (2. * speedOfSoundSq) * (L5 + L1) + timeStep / (2. * speedOfSoundSq) * (T5 + T1);  //!use density according to relaxation approach to enable an exit of acoustic waves
	veloX = veloBoundaryNode[0] - timeStep / (2. * densityBoundaryNode * speedOfSound) * (L5 - L1) + timeStep / (2. * densityBoundaryNode * speedOfSound) * (T5 - T1);
	veloY = veloBoundaryNode[1] - timeStep * L3 + timeStep * T3;
	veloZ = veloBoundaryNode[2] - timeStep * L4 + timeStep * T4;

	velocity[0] = veloX;
	velocity[1] = veloY;
	velocity[2] = veloZ;

}//!end LODIpressureCalcVelocity


///////////////////////////////////////////////////////////////////
//!					SRT	COLLISION OPERATORS
///////////////////////////////////////////////////////////////////

void Node::collideSRT(Control& bc, bool& alternating)
{
	//!calculate density
	double dens = this->calcDensity(bc, alternating);

	//!calculate velocity
	double velo[3];
	this->calcVelocity(bc, alternating, dens, velo);

	//!calculate equilibrium from density and velocity
	double equil[19];

    if (bc.getOrderOfEquilibrium() == 2)
	{
        this->calcEquilibrium2_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 3)
	{
        this->calcEquilibrium3_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 4)
	{
        this->calcEquilibrium4_GHbasis(equil, bc, dens, velo);
	}
	else
	{
        cerr << "Wrong order of equilibrium!" << flush; exit(1);
	}

	for (int dir = 0; dir < 19; dir++)
	{
		distr_[alternating][dir] += colFreq_ * (equil[dir] - distr_[alternating][dir]);
	}

	//!set macroscopic density and velocity
	double densityDummy = this->calcDensity(bc, alternating);
	double veloDummy[3] = { 0. };
	this->calcVelocity(bc, alternating, densityDummy, veloDummy);

	this->setNodeDensity(densityDummy);
	int dirX = 0;
	int dirY = 1;
	int dirZ = 2;
	this->setNodeVelocityComponent(dirX, veloDummy[0]);
	this->setNodeVelocityComponent(dirY, veloDummy[1]);
	this->setNodeVelocityComponent(dirZ, veloDummy[2]);
}


void Node::collideRR(Control& bc, bool& alternating)
{
	//!calculate density
	double dens = this->calcDensity(bc, alternating);

	//!calculate velocity
	double velo[3];
	this->calcVelocity(bc, alternating, dens, velo);

	//!calculate equilibrium from density and velocity
	double equil[19];

    if (bc.getOrderOfEquilibrium() == 2)
	{
        this->calcEquilibrium2_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 3)
	{
        this->calcEquilibrium3_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 4)
	{
        this->calcEquilibrium4_GHbasis(equil, bc, dens, velo);
	}
	else
	{
        cerr << "Wrong order of equilibrium!" << flush; exit(1);
	}

	//!calculate strain-rate tensor through second order central difference scheme (strain-rate tensor is a member of node class)
	//!calculate also cubic Mach correction terms to overcome low symmetry limitations of isothermal lattice (psi is a member of node class)
	if (bc.getCubicMachCorrection() == "yes")
	{
		this->calcStrainRateTensorFDandCubicMachCorrection(bc);
	}
	else
	{

	}

	//!calculate deviatoric stress tensor (deviatoric stress tensor is a member of node class)
	this->calcPi1Tensor(alternating, equil, bc, velo);

	//!calculate first-order off-equilibrium distribution from recursive-regularized approach (third-order hermite polynomial basis)
	double nEquilRR[19];

    if (bc.getOrderOfEquilibrium() == 2)
	{
        this->calcRRoffEquilibrium2_GHbasis(alternating, nEquilRR, equil, bc, velo); //!PR approach
	}
	else if (bc.getOrderOfEquilibrium() == 3)
	{
        this->calcRRoffEquilibrium3_GHbasis(alternating, nEquilRR, equil, bc, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 4)
	{
        this->calcRRoffEquilibrium4_GHbasis(alternating, nEquilRR, equil, bc, velo);
	}
	else
	{
        cerr << "Wrong order of equilibrium!" << flush; exit(1);
	}

	if (bc.getCubicMachCorrection() == "yes")
	{
		for (int dir = 0; dir < 19; dir++)
		{
			distr_[alternating][dir] = equil[dir] + (1 - colFreq_) * nEquilRR[dir] + 0.5 * timeStep_ * psi_[dir];
		}
	}
	else
	{
		for (int dir = 0; dir < 19; dir++)
		{
			distr_[alternating][dir] = equil[dir] + (1 - colFreq_) * nEquilRR[dir];
		}
	}

	//!set macroscopic density and velocity
	double densityDummy = this->calcDensity(bc, alternating);
	double veloDummy[3] = { 0. };
	this->calcVelocity(bc, alternating, densityDummy, veloDummy);

	this->setNodeDensity(densityDummy);
	int dirX = 0;
	int dirY = 1;
	int dirZ = 2;
	this->setNodeVelocityComponent(dirX, veloDummy[0]);
	this->setNodeVelocityComponent(dirY, veloDummy[1]);
	this->setNodeVelocityComponent(dirZ, veloDummy[2]);
}


void Node::collideHRR(Control& bc, bool& alternating)
{
	//!calculate density
	double dens = this->calcDensity(bc, alternating);

	//!calculate velocity
	double velo[3];
	this->calcVelocity(bc, alternating, dens, velo);

	//!calculate equilibrium from density and velocity
	double equil[19];

    if (bc.getOrderOfEquilibrium() == 2)
	{
        this->calcEquilibrium2_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 3)
	{
        this->calcEquilibrium3_GHbasis(equil, bc, dens, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 4)
	{
        this->calcEquilibrium4_GHbasis(equil, bc, dens, velo);
	}
	else
	{
        cerr << "Wrong order of equilibrium!" << flush; exit(1);
	}

	//!calculate strain-rate tensor through second order central difference scheme (strain-rate tensor is a member of node class)
	//!calculate also cubic Mach correction terms to overcome low symmetry limitations of isothermal lattice (psi is a member of node class)
	this->calcStrainRateTensorFDandCubicMachCorrection(bc);

	//!calculate deviatoric stress tensor (deviatoric stress tensor is a member of node class)
	this->calcPi1Tensor(alternating, equil, bc, velo);

	//!calculate first-order off-equilibrium distribution from recursive-regularized approach (third-order hermite polynomial basis)
	double nEquilHRR[19];

    if (bc.getOrderOfEquilibrium() == 2)
	{
        this->calcHRRoffEquilibrium2_GHbasis(alternating, nEquilHRR, equil, bc, velo); //!PR approach
	}
	else if (bc.getOrderOfEquilibrium() == 3)
	{
        this->calcHRRoffEquilibrium3_GHbasis(alternating, nEquilHRR, equil, bc, velo);
	}
	else if (bc.getOrderOfEquilibrium() == 4)
	{
        this->calcHRRoffEquilibrium4_GHbasis(alternating, nEquilHRR, equil, bc, velo);
	}
	else
	{
        cerr << "Wrong order of equilibrium!" << flush; exit(1);
	}

	if (bc.getCubicMachCorrection() == "yes")
	{
		for (int dir = 0; dir < 19; dir++)
		{
			distr_[alternating][dir] = equil[dir] + (1 - colFreq_) * nEquilHRR[dir] + 0.5 * timeStep_ * psi_[dir];
		}
	}
	else
	{
		for (int dir = 0; dir < 19; dir++)
		{
			distr_[alternating][dir] = equil[dir] + (1 - colFreq_) * nEquilHRR[dir];
		}
	}

	//!set macroscopic density and velocity
	double densityDummy = this->calcDensity(bc, alternating);
	double veloDummy[3] = { 0. };
	this->calcVelocity(bc, alternating, densityDummy, veloDummy);

	this->setNodeDensity(densityDummy);
	int dirX = 0;
	int dirY = 1;
	int dirZ = 2;
	this->setNodeVelocityComponent(dirX, veloDummy[0]);
	this->setNodeVelocityComponent(dirY, veloDummy[1]);
	this->setNodeVelocityComponent(dirZ, veloDummy[2]);
}


void Node::calcEquilibrium2_GHbasis(double* equil, Control& bc, double& dens, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double oneOver2 = 1. / 2.;
	double speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
	double speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, squared
	double reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double xsi = 1. / bc.getMolecularVelocity();
	double xsiSq = pow(xsi, 2.);
	double ux = velo[0];
	double uy = velo[1];
	double uz = velo[2];
	double uxSq = pow(ux, 2);
	double uySq = pow(uy, 2);
	double uzSq = pow(uz, 2);

	//!elements of zeroth-order hermite tensor
	double hermZrth;

	//!elements of first-order hermite tensor
	double hermFrstX, hermFrstY, hermFrstZ;

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!zeroth-order coefficients of equilibrium projection to hermite basis
	double aEqZrth = 1.;

	//!first-order coefficients of equilibrium projection to hermite basis
	double aEqFrstX = velo[0] * aEqZrth;
	double aEqFrstY = velo[1] * aEqZrth;
	double aEqFrstZ = velo[2] * aEqZrth;

	//!second-order coefficients of equilibrium projection to hermite basis
	double aEqSndXX = velo[0] * aEqFrstX;
	double aEqSndYY = velo[1] * aEqFrstY;
	double aEqSndZZ = velo[2] * aEqFrstZ;
	double aEqSndXY = velo[0] * aEqFrstY;
	double aEqSndXZ = velo[0] * aEqFrstZ;
	double aEqSndYZ = velo[1] * aEqFrstZ;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate zeroth-order hermite tensor elements
		hermZrth = 1.;

		//!calculate first-order hermite tensor elements
		hermFrstX = xsiX[dir];
		hermFrstY = xsiY[dir];
		hermFrstZ = xsiZ[dir];

		//!calculate second-order hermite tensor elements
		hermSndXX = xsiX[dir] * xsiX[dir] - speedOfSoundSq;
		hermSndYY = xsiY[dir] * xsiY[dir] - speedOfSoundSq;
		hermSndZZ = xsiZ[dir] * xsiZ[dir] - speedOfSoundSq;
		hermSndXY = xsiX[dir] * xsiY[dir];
		hermSndXZ = xsiX[dir] * xsiZ[dir];
		hermSndYZ = xsiY[dir] * xsiZ[dir];

		equil[dir] = dens * weight[dir] * (hermZrth * aEqZrth //!zeroth-
			+ reciSpeedOfSoundSq * (hermFrstX * aEqFrstX //!first-
				+ hermFrstY * aEqFrstY
				+ hermFrstZ * aEqFrstZ)
			+ oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aEqSndXX //!second-order terms from hermite polynomial approximation
				+ hermSndYY * aEqSndYY
				+ hermSndZZ * aEqSndZZ
				+ hermSndXY * aEqSndXY * 2.
				+ hermSndXZ * aEqSndXZ * 2.
				+ hermSndYZ * aEqSndYZ * 2.));
	}
}


void Node::calcEquilibrium3_GHbasis(double* equil, Control& bc, double& dens, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double oneOver2 = 1. / 2.;
	double oneOver6 = 1. / 6.;
	double speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
	double speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double speedOfSoundHc = pow(bc.getSpeedOfSound(), 6.); //!speed of sound, hexic
	double reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, squared
	double reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double reciSpeedOfSoundHc = 1. / speedOfSoundHc; //!reciprocal value of speed of sound, hexic
	double xsi = 1. / bc.getMolecularVelocity();
	double xsiSq = pow(xsi, 2.);
	double xsiPow3 = pow(xsi, 3.);
	double xsiPow4 = pow(xsi, 4.);
	double ux = velo[0];
	double uy = velo[1];
	double uz = velo[2];
	double uxSq = pow(ux, 2.);
	double uySq = pow(uy, 2.);
	double uzSq = pow(uz, 2.);

	//!elements of zeroth-order hermite tensor
	double hermZrth;

	//!elements of first-order hermite tensor
	double hermFrstX, hermFrstY, hermFrstZ;

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!elements of third-order hermite tensor (complete projection not possible due to insufficient quadrature using D3Q19)
	double hermThrdXXY, hermThrdXYY, hermThrdXXZ, hermThrdXZZ, hermThrdYYZ, hermThrdYZZ, hermThrdXXYpYZZ, hermThrdXXZpYYZ, hermThrdXYYpXZZ, hermThrdYZZmXXY, hermThrdYYZmXXZ, hermThrdXZZmXYY;

	//!zeroth-order coefficients of equilibrium projection to hermite basis
	double aEqZrth = 1.;

	//!first-order coefficients of equilibrium projection to hermite basis
	double aEqFrstX = velo[0] * aEqZrth;
	double aEqFrstY = velo[1] * aEqZrth;
	double aEqFrstZ = velo[2] * aEqZrth;

	//!second-order coefficients of equilibrium projection to hermite basis
	double aEqSndXX = aEqFrstX * aEqFrstX;
	double aEqSndYY = aEqFrstY * aEqFrstY;
	double aEqSndZZ = aEqFrstZ * aEqFrstZ;
	double aEqSndXY = aEqFrstX * aEqFrstY;
	double aEqSndXZ = aEqFrstX * aEqFrstZ;
	double aEqSndYZ = aEqFrstY * aEqFrstZ;

	//!third-order coefficients of equilibrium projection hermite basis
	double aEqThrdXXY = aEqFrstY * aEqSndXX;
	double aEqThrdXYY = aEqFrstX * aEqSndYY;
	double aEqThrdXXZ = aEqFrstZ * aEqSndXX;
	double aEqThrdXZZ = aEqFrstX * aEqSndZZ;
	double aEqThrdYYZ = aEqFrstZ * aEqSndYY;
	double aEqThrdYZZ = aEqFrstY * aEqSndZZ;

	//!order 3 fully orthogonal
	double aEqThrdXXYpYZZ = aEqThrdXXY + aEqThrdYZZ;
	double aEqThrdXXZpYYZ = aEqThrdXXZ + aEqThrdYYZ;
	double aEqThrdXYYpXZZ = aEqThrdXYY + aEqThrdXZZ;
	double aEqThrdYZZmXXY = aEqThrdYZZ - aEqThrdXXY;
	double aEqThrdYYZmXXZ = aEqThrdYYZ - aEqThrdXXZ;
	double aEqThrdXZZmXYY = aEqThrdXZZ - aEqThrdXYY;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate zeroth-order hermite tensor elements
		hermZrth = 1.;

		//!calculate first-order hermite tensor elements
		hermFrstX = xsiX[dir];
		hermFrstY = xsiY[dir];
		hermFrstZ = xsiZ[dir];

		//!calculate second-order hermite tensor elements
		hermSndXX = hermFrstX * hermFrstX - speedOfSoundSq;
		hermSndYY = hermFrstY * hermFrstY - speedOfSoundSq;
		hermSndZZ = hermFrstZ * hermFrstZ - speedOfSoundSq;
		hermSndXY = hermFrstX * hermFrstY;
		hermSndXZ = hermFrstX * hermFrstZ;
		hermSndYZ = hermFrstY * hermFrstZ;

		//!calculate third-order hermite tensor elements
		hermThrdXXY = hermSndXX * hermFrstY;
		hermThrdXYY = hermSndYY * hermFrstX;
		hermThrdXXZ = hermSndXX * hermFrstZ;
		hermThrdXZZ = hermSndZZ * hermFrstX;
		hermThrdYYZ = hermSndYY * hermFrstZ;
		hermThrdYZZ = hermSndZZ * hermFrstY;

		//!order 3 fully orthogonal
		hermThrdXXYpYZZ = hermThrdXXY + hermThrdYZZ;
		hermThrdXXZpYYZ = hermThrdXXZ + hermThrdYYZ;
		hermThrdXYYpXZZ = hermThrdXYY + hermThrdXZZ;
		hermThrdYZZmXXY = hermThrdYZZ - hermThrdXXY;
		hermThrdYYZmXXZ = hermThrdYYZ - hermThrdXXZ;
		hermThrdXZZmXYY = hermThrdXZZ - hermThrdXYY;

		equil[dir] = dens * weight[dir] * (hermZrth * aEqZrth //!zeroth-
			+ reciSpeedOfSoundSq * (hermFrstX * aEqFrstX //!first-
				+ hermFrstY * aEqFrstY
				+ hermFrstZ * aEqFrstZ)
			+ oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aEqSndXX //!second-order terms from hermite polynomial approximation
				+ hermSndYY * aEqSndYY
				+ hermSndZZ * aEqSndZZ
				+ hermSndXY * aEqSndXY * 2.
				+ hermSndXZ * aEqSndXZ * 2.
				+ hermSndYZ * aEqSndYZ * 2.)

			//+ oneOver6 * reciSpeedOfSoundHc * (hermThrdXXY * aEqThrdXXY * 3. //!third-order terms from hermite polynomial approximation
			//								 + hermThrdXXZ * aEqThrdXXZ * 3.
			//								 + hermThrdXYY * aEqThrdXYY * 3.
			//								 + hermThrdXZZ * aEqThrdXZZ * 3.
			//								 + hermThrdYYZ * aEqThrdYYZ * 3.
			//								 + hermThrdYZZ * aEqThrdYZZ * 3.)

			//!to avoid spurious couplings and thus stability issues for the D3Q19-GH formalism rely on fully orthogonal basis
			//!orthogonalize fourth-order hermite polynomials with second-order ones and further orthogonalize them between each other
			//!see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. G17 - G19

			//!orthogonal third-order hermite polynomials (see Jacob et al. 2018, Journal of Turbulence)
			+oneOver6 * reciSpeedOfSoundHc * (hermThrdXXYpYZZ * aEqThrdXXYpYZZ * 3.
				+ hermThrdXXZpYYZ * aEqThrdXXZpYYZ * 3.
				+ hermThrdXYYpXZZ * aEqThrdXYYpXZZ * 3.
				+ hermThrdYZZmXXY * aEqThrdYZZmXXY
				+ hermThrdYYZmXXZ * aEqThrdYYZmXXZ
				+ hermThrdXZZmXYY * aEqThrdXZZmXYY));
	}
}


void Node::calcEquilibrium4_GHbasis(double* equil, Control& bc, double& dens, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double oneOver2 = 1. / 2.;
	double oneOver4 = 1. / 4.;
	double oneOver6 = 1. / 6.;
	double twoOver7 = 2. / 7.;
	double twoOver5 = 2. / 5.;
	double eightOver7 = 8. / 7.;
	double fiftysixOver45 = 56. / 45.;
	double fortyOver27 = 40. / 27.;
	double speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
	double speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double speedOfSoundHc = pow(bc.getSpeedOfSound(), 6.); //!speed of sound, hexic
	double speedOfSoundOc = pow(bc.getSpeedOfSound(), 8.); //!speed of sound, octic
	double reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, squared
	double reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double reciSpeedOfSoundHc = 1. / speedOfSoundHc; //!reciprocal value of speed of sound, hexic
	double reciSpeedOfSoundOc = 1. / speedOfSoundOc; //!reciprocal value of speed of sound, octic
	double xsi = 1. / bc.getMolecularVelocity();
	double xsiSq = pow(xsi, 2.);
	double xsiPow3 = pow(xsi, 3.);
	double xsiPow4 = pow(xsi, 4.);
	double ux = velo[0];
	double uy = velo[1];
	double uz = velo[2];
	double uxSq = pow(ux, 2.);
	double uySq = pow(uy, 2.);
	double uzSq = pow(uz, 2.);

	//!elements of zeroth-order hermite tensor
	double hermZrth;

	//!elements of first-order hermite tensor
	double hermFrstX, hermFrstY, hermFrstZ;

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!elements of third-order hermite tensor (complete projection not possible due to insufficient quadrature using D3Q19)
	double hermThrdXXY, hermThrdXYY, hermThrdXXZ, hermThrdXZZ, hermThrdYYZ, hermThrdYZZ, hermThrdXXYpYZZ, hermThrdXXZpYYZ, hermThrdXYYpXZZ, hermThrdYZZmXXY, hermThrdYYZmXXZ, hermThrdXZZmXYY;

	//!elements of fourth-order hermite tensor (complete projection not possible due to insufficient quadrature using D3Q19)
	double hermFrthXXYY, hermFrthXXZZ, hermFrthYYZZ, hermFrthXXYY4o2, hermFrthXXZZ4o2, hermFrthYYZZ4o2, hermFrthXXYYFO, hermFrthXXZZFO, hermFrthYYZZFO;

	//!zeroth-order coefficients of equilibrium projection to hermite basis
	double aEqZrth = 1.;

	//!first-order coefficients of equilibrium projection to hermite basis
	double aEqFrstX = velo[0] * aEqZrth;
	double aEqFrstY = velo[1] * aEqZrth;
	double aEqFrstZ = velo[2] * aEqZrth;

	//!second-order coefficients of equilibrium projection to hermite basis
	double aEqSndXX = aEqFrstX * aEqFrstX;
	double aEqSndYY = aEqFrstY * aEqFrstY;
	double aEqSndZZ = aEqFrstZ * aEqFrstZ;
	double aEqSndXY = aEqFrstX * aEqFrstY;
	double aEqSndXZ = aEqFrstX * aEqFrstZ;
	double aEqSndYZ = aEqFrstY * aEqFrstZ;

	//!third-order coefficients of equilibrium projection hermite basis
	double aEqThrdXXY = aEqFrstY * aEqSndXX;
	double aEqThrdXYY = aEqFrstX * aEqSndYY;
	double aEqThrdXXZ = aEqFrstZ * aEqSndXX;
	double aEqThrdXZZ = aEqFrstX * aEqSndZZ;
	double aEqThrdYYZ = aEqFrstZ * aEqSndYY;
	double aEqThrdYZZ = aEqFrstY * aEqSndZZ;

	//!order 3 fully orthogonal
	double aEqThrdXXYpYZZ = aEqThrdXXY + aEqThrdYZZ;
	double aEqThrdXXZpYYZ = aEqThrdXXZ + aEqThrdYYZ;
	double aEqThrdXYYpXZZ = aEqThrdXYY + aEqThrdXZZ;
	double aEqThrdYZZmXXY = aEqThrdYZZ - aEqThrdXXY;
	double aEqThrdYYZmXXZ = aEqThrdYYZ - aEqThrdXXZ;
	double aEqThrdXZZmXYY = aEqThrdXZZ - aEqThrdXYY;

	//!fourth-order coefficients of equilibrium projection to hermite basis
	double aEqFrthXXYY = aEqFrstY * aEqThrdXXY;
	double aEqFrthXXZZ = aEqFrstZ * aEqThrdXXZ;
	double aEqFrthYYZZ = aEqFrstY * aEqThrdYZZ;

	//!order 4 orthogonal with respect to order 2
	double aEqFrthXXYY4o2 = aEqFrthXXYY + oneOver6 * aEqSndZZ;
	double aEqFrthXXZZ4o2 = aEqFrthXXZZ + oneOver6 * aEqSndYY;
	double aEqFrthYYZZ4o2 = aEqFrthYYZZ + oneOver6 * aEqSndXX;

	//!order 4 fully orthogonal
	double aEqFrthXXYYFO = aEqFrthXXYY4o2;
	double aEqFrthXXZZFO = aEqFrthXXZZ4o2 + twoOver7 * aEqFrthXXYYFO;
	double aEqFrthYYZZFO = aEqFrthYYZZ4o2 + twoOver7 * aEqFrthXXYYFO + twoOver5 * aEqFrthXXZZFO;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate zeroth-order hermite tensor elements
		hermZrth = 1.;

		//!calculate first-order hermite tensor elements
		hermFrstX = xsiX[dir];
		hermFrstY = xsiY[dir];
		hermFrstZ = xsiZ[dir];

		//!calculate second-order hermite tensor elements
		hermSndXX = hermFrstX * hermFrstX - speedOfSoundSq;
		hermSndYY = hermFrstY * hermFrstY - speedOfSoundSq;
		hermSndZZ = hermFrstZ * hermFrstZ - speedOfSoundSq;
		hermSndXY = hermFrstX * hermFrstY;
		hermSndXZ = hermFrstX * hermFrstZ;
		hermSndYZ = hermFrstY * hermFrstZ;

		//!calculate third-order hermite tensor elements
		hermThrdXXY = hermSndXX * hermFrstY;
		hermThrdXYY = hermSndYY * hermFrstX;
		hermThrdXXZ = hermSndXX * hermFrstZ;
		hermThrdXZZ = hermSndZZ * hermFrstX;
		hermThrdYYZ = hermSndYY * hermFrstZ;
		hermThrdYZZ = hermSndZZ * hermFrstY;

		//!order 3 fully orthogonal
		hermThrdXXYpYZZ = hermThrdXXY + hermThrdYZZ;
		hermThrdXXZpYYZ = hermThrdXXZ + hermThrdYYZ;
		hermThrdXYYpXZZ = hermThrdXYY + hermThrdXZZ;
		hermThrdYZZmXXY = hermThrdYZZ - hermThrdXXY;
		hermThrdYYZmXXZ = hermThrdYYZ - hermThrdXXZ;
		hermThrdXZZmXYY = hermThrdXZZ - hermThrdXYY;

		//!calculate fourth-order hermite tensor elements
		hermFrthXXYY = hermSndXX * hermSndYY;
		hermFrthXXZZ = hermSndXX * hermSndZZ;
		hermFrthYYZZ = hermSndYY * hermSndZZ;

		//!order 4 orthogonal with respect to order 2
		hermFrthXXYY4o2 = hermFrthXXYY + oneOver6 * hermSndZZ;
		hermFrthXXZZ4o2 = hermFrthXXZZ + oneOver6 * hermSndYY;
		hermFrthYYZZ4o2 = hermFrthYYZZ + oneOver6 * hermSndXX;

		//!order 4 fully orthogonal
		hermFrthXXYYFO = hermFrthXXYY4o2;
		hermFrthXXZZFO = hermFrthXXZZ4o2 + twoOver7 * hermFrthXXYYFO;
		hermFrthYYZZFO = hermFrthYYZZ4o2 + twoOver7 * hermFrthXXYYFO + twoOver5 * hermFrthXXZZFO;

		equil[dir] = dens * weight[dir] * (hermZrth * aEqZrth //!zeroth-
			+ reciSpeedOfSoundSq * (hermFrstX * aEqFrstX //!first-
				+ hermFrstY * aEqFrstY
				+ hermFrstZ * aEqFrstZ)
			+ oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aEqSndXX //!second-order terms from hermite polynomial approximation
				+ hermSndYY * aEqSndYY
				+ hermSndZZ * aEqSndZZ
				+ hermSndXY * aEqSndXY * 2.
				+ hermSndXZ * aEqSndXZ * 2.
				+ hermSndYZ * aEqSndYZ * 2.)

			//+ oneOver6 * reciSpeedOfSoundHc * (hermThrdXXY * aEqThrdXXY * 3. //!third-order terms from hermite polynomial approximation
			//								 + hermThrdXXZ * aEqThrdXXZ * 3.
			//								 + hermThrdXYY * aEqThrdXYY * 3.
			//								 + hermThrdXZZ * aEqThrdXZZ * 3.
			//								 + hermThrdYYZ * aEqThrdYYZ * 3.
			//								 + hermThrdYZZ * aEqThrdYZZ * 3.)
			//+ oneOver4 * reciSpeedOfSoundOc * (hermFrthXXYY * aEqFrthXXYY //!fourth-order terms from hermite polynomial approximation
			//							     + hermFrthXXZZ * aEqFrthXXZZ
			//							     + hermFrthYYZZ * aEqFrthYYZZ));

			//!to avoid spurious couplings and thus stability issues for the D3Q19-GH formalism rely on fully orthogonal basis
			//!orthogonalize fourth-order hermite polynomials with second-order ones and further orthogonalize them between each other
			//!see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. G17 - G19

			//!orthogonal third-order hermite polynomials (see Jacob et al. 2018, Journal of Turbulence)
			+oneOver6 * reciSpeedOfSoundHc * (hermThrdXXYpYZZ * aEqThrdXXYpYZZ * 3.
				+ hermThrdXXZpYYZ * aEqThrdXXZpYYZ * 3.
				+ hermThrdXYYpXZZ * aEqThrdXYYpXZZ * 3.
				+ hermThrdYZZmXXY * aEqThrdYZZmXXY
				+ hermThrdYYZmXXZ * aEqThrdYYZmXXZ
				+ hermThrdXZZmXYY * aEqThrdXZZmXYY)

		//!fully orthogonal fourth-order hermite polynomials
	   + oneOver4 * reciSpeedOfSoundOc * (eightOver7 * hermFrthXXYYFO * aEqFrthXXYYFO + fiftysixOver45 * hermFrthXXZZFO * aEqFrthXXZZFO + fortyOver27 * hermFrthYYZZFO * aEqFrthYYZZFO));
	}
}



void Node::calcRRoffEquilibrium2_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
    double  speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double  reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, quartic
    double	reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double  oneOver2 = 1. / 2.;

	//!ATTENTION: FIRST ORDER TERMS OF NON-EQUILIBRIUM PROJECTION TO HERMITE-BASIS ARE NEEDED WITH GUO FORCING! (SEE: 10.1016/J.compfluid.2019.06.030)
	//!elements of first-order hermite tensor
	double hermFrstX[19] = { 0. }, hermFrstY[19] = { 0. }, hermFrstZ[19] = { 0. };

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!coefficients of first-order non-equilibrium projection to second-order hermite basis
	double aNeqSndXX = Pi1Tensor_[0][0];
	double aNeqSndYY = Pi1Tensor_[1][1];
	double aNeqSndZZ = Pi1Tensor_[2][2];
	double aNeqSndXY = Pi1Tensor_[0][1];
	double aNeqSndXZ = Pi1Tensor_[0][2];
	double aNeqSndYZ = Pi1Tensor_[1][2];

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate second-order hermite tensor elements
		hermSndXX = xsiX[dir] * xsiX[dir] - speedOfSoundSq;
		hermSndYY = xsiY[dir] * xsiY[dir] - speedOfSoundSq;
		hermSndZZ = xsiZ[dir] * xsiZ[dir] - speedOfSoundSq;
		hermSndXY = xsiX[dir] * xsiY[dir];
		hermSndXZ = xsiX[dir] * xsiZ[dir];
		hermSndYZ = xsiY[dir] * xsiZ[dir];

		//reciSpeedOfSoundSq * (hermFrstX[dir] * aNeqFrstX + hermFrstY[dir] * aNeqFrstY + hermFrstZ[dir] * aNeqFrstZ) +

		nEquilRR[dir] = weight[dir] * (oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aNeqSndXX + hermSndYY * aNeqSndYY + hermSndZZ * aNeqSndZZ
			+ hermSndXY * aNeqSndXY * 2. + hermSndXZ * aNeqSndXZ * 2. + hermSndYZ * aNeqSndYZ * 2.)); //!second-order terms from hermite polynomial approximation
	}
}




void Node::calcRRoffEquilibrium3_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
	double  speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double  speedOfSoundHc = pow(bc.getSpeedOfSound(), 6.); //!speed of sound, hexic
	double  reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundHc = 1. / speedOfSoundHc; //!reciprocal value of speed of sound, hexic
	double  oneOver2 = 1. / 2.;
	double  oneOver6 = 1. / 6.;

	//!ATTENTION: FIRST ORDER TERMS OF NON-EQUILIBRIUM PROJECTION TO HERMITE-BASIS ARE NEEDED WITH GUO FORCING! (SEE: 10.1016/J.compfluid.2019.06.030)
	//!elements of first-order hermite tensor
	double hermFrstX[19] = { 0. }, hermFrstY[19] = { 0. }, hermFrstZ[19] = { 0. };

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!elements of third-order hermite tensor
	double hermThrdXXY, hermThrdXYY, hermThrdXXZ, hermThrdXZZ, hermThrdYYZ, hermThrdYZZ, hermThrdXXYpYZZ, hermThrdXXZpYYZ, hermThrdXYYpXZZ, hermThrdYZZmXXY, hermThrdYYZmXXZ, hermThrdXZZmXYY;

	//!coefficients of first-order non-equilibrium projection to second-order hermite basis
	double aNeqSndXX = Pi1Tensor_[0][0];
	double aNeqSndYY = Pi1Tensor_[1][1];
	double aNeqSndZZ = Pi1Tensor_[2][2];
	double aNeqSndXY = Pi1Tensor_[0][1];
	double aNeqSndXZ = Pi1Tensor_[0][2];
	double aNeqSndYZ = Pi1Tensor_[1][2];

	//!coefficients of first-order non-equilibrium projection to third-order hermite basis via their recursive properties
	double aNeqThrdXXY = 2. * velo[0] * aNeqSndXY + velo[1] * aNeqSndXX;
	double aNeqThrdXYY = 2. * velo[1] * aNeqSndXY + velo[0] * aNeqSndYY;
	double aNeqThrdXXZ = 2. * velo[0] * aNeqSndXZ + velo[2] * aNeqSndXX;
	double aNeqThrdXZZ = 2. * velo[2] * aNeqSndXZ + velo[0] * aNeqSndZZ;
	double aNeqThrdYYZ = 2. * velo[1] * aNeqSndYZ + velo[2] * aNeqSndYY;
	double aNeqThrdYZZ = 2. * velo[2] * aNeqSndYZ + velo[1] * aNeqSndZZ;

	//!orthogonalize third-order hermite polynomials between each other
	double aNeqThrdXXYpYZZ = aNeqThrdXXY + aNeqThrdYZZ;
	double aNeqThrdXXZpYYZ = aNeqThrdXXZ + aNeqThrdYYZ;
	double aNeqThrdXYYpXZZ = aNeqThrdXYY + aNeqThrdXZZ;
	double aNeqThrdYZZmXXY = aNeqThrdYZZ - aNeqThrdXXY;
	double aNeqThrdYYZmXXZ = aNeqThrdYYZ - aNeqThrdXXZ;
	double aNeqThrdXZZmXYY = aNeqThrdXZZ - aNeqThrdXYY;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate second-order hermite tensor elements
		hermSndXX = xsiX[dir] * xsiX[dir] - speedOfSoundSq;
		hermSndYY = xsiY[dir] * xsiY[dir] - speedOfSoundSq;
		hermSndZZ = xsiZ[dir] * xsiZ[dir] - speedOfSoundSq;
		hermSndXY = xsiX[dir] * xsiY[dir];
		hermSndXZ = xsiX[dir] * xsiZ[dir];
		hermSndYZ = xsiY[dir] * xsiZ[dir];

		//!calculate third-order hermite tensor elements
		hermThrdXXY = hermSndXX * xsiY[dir];
		hermThrdXYY = hermSndYY * xsiX[dir];
		hermThrdXXZ = hermSndXX * xsiZ[dir];
		hermThrdXZZ = hermSndZZ * xsiX[dir];
		hermThrdYYZ = hermSndYY * xsiZ[dir];
		hermThrdYZZ = hermSndZZ * xsiY[dir];

		//!orthogonalize third-order hermite polynomials between each other
		hermThrdXXYpYZZ = hermThrdXXY + hermThrdYZZ;
		hermThrdXXZpYYZ = hermThrdXXZ + hermThrdYYZ;
		hermThrdXYYpXZZ = hermThrdXYY + hermThrdXZZ;
		hermThrdYZZmXXY = hermThrdYZZ - hermThrdXXY;
		hermThrdYYZmXXZ = hermThrdYYZ - hermThrdXXZ;
		hermThrdXZZmXYY = hermThrdXZZ - hermThrdXYY;

		nEquilRR[dir] = weight[dir] * (oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aNeqSndXX + hermSndYY * aNeqSndYY + hermSndZZ * aNeqSndZZ
			+ hermSndXY * aNeqSndXY * 2. + hermSndXZ * aNeqSndXZ * 2. + hermSndYZ * aNeqSndYZ * 2.) //!second-order terms from hermite polynomial approximation

//+ oneOver6 * reciSpeedOfSoundHc * (hermThrdXXY * aNeqThrdXXY * 3.
   //							  + hermThrdXYY * aNeqThrdXYY * 3.
   //							  + hermThrdXXZ * aNeqThrdXXZ * 3.
   //							  + hermThrdXZZ * aNeqThrdXZZ * 3.
   //						      + hermThrdYYZ * aNeqThrdYYZ * 3.
   //							  + hermThrdYZZ * aNeqThrdYZZ * 3.) //!third-order terms from hermite polynomial approximation

 //!to avoid spurious couplings and thus stability issues for the D3Q19-GH formalism rely on fully orthogonal basis
 //!orthogonalize fourth-order hermite polynomials with second-order ones and further orthogonalize them between each other
 //see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. G17 - G19

 //!orthogonal third-order hermite polynomials (see Jacob et al. 2018, Journal of Turbulence)
			+oneOver6 * reciSpeedOfSoundHc * (hermThrdXXYpYZZ * aNeqThrdXXYpYZZ * 3.
				+ hermThrdXXZpYYZ * aNeqThrdXXZpYYZ * 3.
				+ hermThrdXYYpXZZ * aNeqThrdXYYpXZZ * 3.
				+ hermThrdYZZmXXY * aNeqThrdYZZmXXY
				+ hermThrdYYZmXXZ * aNeqThrdYYZmXXZ
				+ hermThrdXZZmXYY * aNeqThrdXZZmXYY));
	}
}



void Node::calcRRoffEquilibrium4_GHbasis(bool& alternating, double* nEquilRR, double* equil, Control& bc, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
	double  speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double  speedOfSoundHc = pow(bc.getSpeedOfSound(), 6.); //!speed of sound, hexic
	double  speedOfSoundOc = pow(bc.getSpeedOfSound(), 8.); //!speed of sound, octic
	double  reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundHc = 1. / speedOfSoundHc; //!reciprocal value of speed of sound, hexic
	double  reciSpeedOfSoundOc = 1. / speedOfSoundOc; //!reciprocal value of speed of sound, octic
	double  oneOver2 = 1. / 2.;
	double  oneOver4 = 1. / 4.;
	double  oneOver6 = 1. / 6.;
	double  twoOver7 = 2. / 7.;
	double  twoOver5 = 2. / 5.;
	double  eightOver7 = 8. / 7.;
	double  fiftysixOver45 = 56. / 45.;
	double  fortyOver27 = 40. / 27.;

	//!ATTENTION: FIRST ORDER TERMS OF NON-EQUILIBRIUM PROJECTION TO HERMITE-BASIS ARE NEEDED WITH GUO FORCING! (SEE: 10.1016/J.compfluid.2019.06.030)
	//!elements of first-order hermite tensor
	double hermFrstX[19] = { 0. }, hermFrstY[19] = { 0. }, hermFrstZ[19] = { 0. };

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!elements of third-order hermite tensor
	double hermThrdXXY, hermThrdXYY, hermThrdXXZ, hermThrdXZZ, hermThrdYYZ, hermThrdYZZ, hermThrdXXYpYZZ, hermThrdXXZpYYZ, hermThrdXYYpXZZ, hermThrdYZZmXXY, hermThrdYYZmXXZ, hermThrdXZZmXYY;

	//!elements of fourth-order hermite tensor (complete projection not possible due to insufficient quadrature using D3Q19)
	double hermFrthXXYY, hermFrthXXZZ, hermFrthYYZZ, hermFrthXXYY4o2, hermFrthXXZZ4o2, hermFrthYYZZ4o2, hermFrthXXYYFO, hermFrthXXZZFO, hermFrthYYZZFO;

	//!coefficients of first-order non-equilibrium projection to second-order hermite basis
	double aNeqSndXX = Pi1Tensor_[0][0];
	double aNeqSndYY = Pi1Tensor_[1][1];
	double aNeqSndZZ = Pi1Tensor_[2][2];
	double aNeqSndXY = Pi1Tensor_[0][1];
	double aNeqSndXZ = Pi1Tensor_[0][2];
	double aNeqSndYZ = Pi1Tensor_[1][2];

	//!coefficients of first-order non-equilibrium projection to third-order hermite basis via their recursive properties
	double aNeqThrdXXY = 2. * velo[0] * aNeqSndXY + velo[1] * aNeqSndXX;
	double aNeqThrdXYY = 2. * velo[1] * aNeqSndXY + velo[0] * aNeqSndYY;
	double aNeqThrdXXZ = 2. * velo[0] * aNeqSndXZ + velo[2] * aNeqSndXX;
	double aNeqThrdXZZ = 2. * velo[2] * aNeqSndXZ + velo[0] * aNeqSndZZ;
	double aNeqThrdYYZ = 2. * velo[1] * aNeqSndYZ + velo[2] * aNeqSndYY;
	double aNeqThrdYZZ = 2. * velo[2] * aNeqSndYZ + velo[1] * aNeqSndZZ;

	//!orthogonalize third-order hermite polynomials between each other
	double aNeqThrdXXYpYZZ = aNeqThrdXXY + aNeqThrdYZZ;
	double aNeqThrdXXZpYYZ = aNeqThrdXXZ + aNeqThrdYYZ;
	double aNeqThrdXYYpXZZ = aNeqThrdXYY + aNeqThrdXZZ;
	double aNeqThrdYZZmXXY = aNeqThrdYZZ - aNeqThrdXXY;
	double aNeqThrdYYZmXXZ = aNeqThrdYYZ - aNeqThrdXXZ;
	double aNeqThrdXZZmXYY = aNeqThrdXZZ - aNeqThrdXYY;

	//!coefficients of first-order non-equilibrium projection to fourth-order hermite basis via their recursive properties
	double aNeqFrthXXYY = velo[1] * velo[1] * aNeqSndXX + velo[0] * velo[0] * aNeqSndYY + 4. * velo[0] * velo[1] * aNeqSndXY;
	double aNeqFrthXXZZ = velo[2] * velo[2] * aNeqSndXX + velo[0] * velo[0] * aNeqSndZZ + 4. * velo[0] * velo[2] * aNeqSndXZ;
	double aNeqFrthYYZZ = velo[1] * velo[1] * aNeqSndZZ + velo[2] * velo[2] * aNeqSndYY + 4. * velo[1] * velo[2] * aNeqSndYZ;

	//!order 4 orthogonal with respect to order 2
	double aNeqFrthXXYY4o2 = aNeqFrthXXYY + oneOver6 * aNeqSndZZ;
	double aNeqFrthXXZZ4o2 = aNeqFrthXXZZ + oneOver6 * aNeqSndYY;
	double aNeqFrthYYZZ4o2 = aNeqFrthYYZZ + oneOver6 * aNeqSndXX;

	//!order 4 fully orthogonal
	double aNeqFrthXXYYFO = aNeqFrthXXYY4o2;
	double aNeqFrthXXZZFO = aNeqFrthXXZZ4o2 + twoOver7 * aNeqFrthXXYYFO;
	double aNeqFrthYYZZFO = aNeqFrthYYZZ4o2 + twoOver7 * aNeqFrthXXYYFO + twoOver5 * aNeqFrthXXZZFO;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate second-order hermite tensor elements
		hermSndXX = xsiX[dir] * xsiX[dir] - speedOfSoundSq;
		hermSndYY = xsiY[dir] * xsiY[dir] - speedOfSoundSq;
		hermSndZZ = xsiZ[dir] * xsiZ[dir] - speedOfSoundSq;
		hermSndXY = xsiX[dir] * xsiY[dir];
		hermSndXZ = xsiX[dir] * xsiZ[dir];
		hermSndYZ = xsiY[dir] * xsiZ[dir];

		//!calculate third-order hermite tensor elements
		hermThrdXXY = hermSndXX * xsiY[dir];
		hermThrdXYY = hermSndYY * xsiX[dir];
		hermThrdXXZ = hermSndXX * xsiZ[dir];
		hermThrdXZZ = hermSndZZ * xsiX[dir];
		hermThrdYYZ = hermSndYY * xsiZ[dir];
		hermThrdYZZ = hermSndZZ * xsiY[dir];

		//!orthogonalize third-order hermite polynomials between each other
		hermThrdXXYpYZZ = hermThrdXXY + hermThrdYZZ;
		hermThrdXXZpYYZ = hermThrdXXZ + hermThrdYYZ;
		hermThrdXYYpXZZ = hermThrdXYY + hermThrdXZZ;
		hermThrdYZZmXXY = hermThrdYZZ - hermThrdXXY;
		hermThrdYYZmXXZ = hermThrdYYZ - hermThrdXXZ;
		hermThrdXZZmXYY = hermThrdXZZ - hermThrdXYY;

		//!calculate fourth-order hermite tensor elements
		hermFrthXXYY = hermSndXX * hermSndYY;
		hermFrthXXZZ = hermSndXX * hermSndZZ;
		hermFrthYYZZ = hermSndYY * hermSndZZ;

		//!order 4 orthogonal with respect to order 2
		hermFrthXXYY4o2 = hermFrthXXYY + oneOver6 * hermSndZZ;
		hermFrthXXZZ4o2 = hermFrthXXZZ + oneOver6 * hermSndYY;
		hermFrthYYZZ4o2 = hermFrthYYZZ + oneOver6 * hermSndXX;

		//!order 4 fully orthogonal
		hermFrthXXYYFO = hermFrthXXYY4o2;
		hermFrthXXZZFO = hermFrthXXZZ4o2 + twoOver7 * hermFrthXXYYFO;
		hermFrthYYZZFO = hermFrthYYZZ4o2 + twoOver7 * hermFrthXXYYFO + twoOver5 * hermFrthXXZZFO;

		nEquilRR[dir] = weight[dir] * (oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aNeqSndXX + hermSndYY * aNeqSndYY + hermSndZZ * aNeqSndZZ
			+ hermSndXY * aNeqSndXY * 2. + hermSndXZ * aNeqSndXZ * 2. + hermSndYZ * aNeqSndYZ * 2.) //!second-order terms from hermite polynomial approximation

//+ oneOver6 * reciSpeedOfSoundHc * (hermThrdXXY * aNeqThrdXXY * 3.
   //							  + hermThrdXYY * aNeqThrdXYY * 3.
   //							  + hermThrdXXZ * aNeqThrdXXZ * 3.
   //							  + hermThrdXZZ * aNeqThrdXZZ * 3.
   //						      + hermThrdYYZ * aNeqThrdYYZ * 3.
   //							  + hermThrdYZZ * aNeqThrdYZZ * 3.) //!third-order terms from hermite polynomial approximation
//+ oneOver4 * reciSpeedOfSoundOc * (hermFrthXXYY * aNeqFrthXXYY + hermFrthXXZZ * aNeqFrthXXZZ + hermFrthYYZZ * aNeqFrthYYZZ)); //!fourth-order terms from hermite polynomial approximation


 //!to avoid spurious couplings and thus stability issues for the D3Q19-GH formalism rely on fully orthogonal basis
 //!orthogonalize fourth-order hermite polynomials with second-order ones and further orthogonalize them between each other
 //see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. G17 - G19

 //!orthogonal third-order hermite polynomials (see Jacob et al. 2018, Journal of Turbulence)
			+oneOver6 * reciSpeedOfSoundHc * (hermThrdXXYpYZZ * aNeqThrdXXYpYZZ * 3.
				+ hermThrdXXZpYYZ * aNeqThrdXXZpYYZ * 3.
				+ hermThrdXYYpXZZ * aNeqThrdXYYpXZZ * 3.
				+ hermThrdYZZmXXY * aNeqThrdYZZmXXY
				+ hermThrdYYZmXXZ * aNeqThrdYYZmXXZ
				+ hermThrdXZZmXYY * aNeqThrdXZZmXYY)

		//!fully orthogonal fourth-order hermite polynomials
		+ oneOver4 * reciSpeedOfSoundOc * (eightOver7 * hermFrthXXYYFO * aNeqFrthXXYYFO + fiftysixOver45 * hermFrthXXZZFO * aNeqFrthXXZZFO + fortyOver27 * hermFrthYYZZFO * aNeqFrthYYZZFO));
	}
}


void Node::calcHRRoffEquilibrium2_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
    double  speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double  reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, quartic
    double	reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double  oneOver2 = 1. / 2.;
	double  hybridPar = bc.getHybridParameter();
	double	timeStep = bc.getTimeStep() / pow(2., level_);

	//!calculate density at node
	double dens = this->calcDensity(bc, alternating);

	//!calculate dynamic viscosity and inverse pre-factor
	double kinVisc = bc.getKinematicViscosity();
	double dynVisc = kinVisc * dens;
	double invPreFactor = -((timeStep * speedOfSoundSq) / (colFreq_ * kinVisc));

	//!ATTENTION: FIRST ORDER TERMS OF NON-EQUILIBRIUM PROJECTION TO HERMITE-BASIS ARE NEEDED WITH GUO FORCING! (SEE: 10.1016/J.compfluid.2019.06.030)
	//!elements of first-order hermite tensor
	double hermFrstX[19] = { 0. }, hermFrstY[19] = { 0. }, hermFrstZ[19] = { 0. };

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!hybrid coefficients including finite-difference approximated strain-rate tensor
	double aHybridXX = hybridPar * Pi1Tensor_[0][0] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][0]);
	double aHybridYY = hybridPar * Pi1Tensor_[1][1] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[1][1]);
	double aHybridZZ = hybridPar * Pi1Tensor_[2][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[2][2]);
	double aHybridXY = hybridPar * Pi1Tensor_[0][1] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][1]);
	double aHybridXZ = hybridPar * Pi1Tensor_[0][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][2]);
	double aHybridYZ = hybridPar * Pi1Tensor_[1][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[1][2]);

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate second-order hermite tensor elements
		hermSndXX = xsiX[dir] * xsiX[dir] - speedOfSoundSq;
		hermSndYY = xsiY[dir] * xsiY[dir] - speedOfSoundSq;
		hermSndZZ = xsiZ[dir] * xsiZ[dir] - speedOfSoundSq;
		hermSndXY = xsiX[dir] * xsiY[dir];
		hermSndXZ = xsiX[dir] * xsiZ[dir];
		hermSndYZ = xsiY[dir] * xsiZ[dir];

		nEquilHRR[dir] = weight[dir] * (oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aHybridXX + hermSndYY * aHybridYY + hermSndZZ * aHybridZZ
			+ hermSndXY * aHybridXY * 2. + hermSndXZ * aHybridXZ * 2. + hermSndYZ * aHybridYZ * 2.)); //!second-order terms from hermite polynomial approximation
	}
}


void Node::calcHRRoffEquilibrium3_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
	double  speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double  speedOfSoundHc = pow(bc.getSpeedOfSound(), 6.); //!speed of sound, hexic
	double  reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundHc = 1. / speedOfSoundHc; //!reciprocal value of speed of sound, hexic
	double  oneOver2 = 1. / 2.;
	double  oneOver6 = 1. / 6.;
	double  hybridPar = bc.getHybridParameter();
	double	timeStep = bc.getTimeStep() / pow(2., level_);

	//!calculate density at node
	double dens = this->calcDensity(bc, alternating);

	//!calculate dynamic viscosity and inverse pre-factor
	double kinVisc = bc.getKinematicViscosity();
	double dynVisc = kinVisc * dens;
	double invPreFactor = -((timeStep * speedOfSoundSq) / (colFreq_ * kinVisc));

	//!ATTENTION: FIRST ORDER TERMS OF NON-EQUILIBRIUM PROJECTION TO HERMITE-BASIS ARE NEEDED WITH GUO FORCING! (SEE: 10.1016/J.compfluid.2019.06.030)
	//!elements of first-order hermite tensor
	double hermFrstX[19] = { 0. }, hermFrstY[19] = { 0. }, hermFrstZ[19] = { 0. };

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!elements of third-order hermite tensor
	double hermThrdXXY, hermThrdXYY, hermThrdXXZ, hermThrdXZZ, hermThrdYYZ, hermThrdYZZ, hermThrdXXYpYZZ, hermThrdXXZpYYZ, hermThrdXYYpXZZ, hermThrdYZZmXXY, hermThrdYYZmXXZ, hermThrdXZZmXYY;

	//!hybrid coefficients including finite-difference approximated strain-rate tensor
	double aHybridXX = hybridPar * Pi1Tensor_[0][0] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][0]);
	double aHybridYY = hybridPar * Pi1Tensor_[1][1] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[1][1]);
	double aHybridZZ = hybridPar * Pi1Tensor_[2][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[2][2]);
	double aHybridXY = hybridPar * Pi1Tensor_[0][1] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][1]);
	double aHybridXZ = hybridPar * Pi1Tensor_[0][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][2]);
	double aHybridYZ = hybridPar * Pi1Tensor_[1][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[1][2]);

	//!coefficients of first-order non-equilibrium projection to third-order hermite basis
	double aNeqThrdXXY = 2. * velo[0] * aHybridXY + velo[1] * aHybridXX;
	double aNeqThrdXYY = 2. * velo[1] * aHybridXY + velo[0] * aHybridYY;
	double aNeqThrdXXZ = 2. * velo[0] * aHybridXZ + velo[2] * aHybridXX;
	double aNeqThrdXZZ = 2. * velo[2] * aHybridXZ + velo[0] * aHybridZZ;
	double aNeqThrdYYZ = 2. * velo[1] * aHybridYZ + velo[2] * aHybridYY;
	double aNeqThrdYZZ = 2. * velo[2] * aHybridYZ + velo[1] * aHybridZZ;

	//!orthogonalize third-order hermite polynomials between each other
	double aNeqThrdXXYpYZZ = aNeqThrdXXY + aNeqThrdYZZ;
	double aNeqThrdXXZpYYZ = aNeqThrdXXZ + aNeqThrdYYZ;
	double aNeqThrdXYYpXZZ = aNeqThrdXYY + aNeqThrdXZZ;
	double aNeqThrdYZZmXXY = aNeqThrdYZZ - aNeqThrdXXY;
	double aNeqThrdYYZmXXZ = aNeqThrdYYZ - aNeqThrdXXZ;
	double aNeqThrdXZZmXYY = aNeqThrdXZZ - aNeqThrdXYY;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate second-order hermite tensor elements
		hermSndXX = xsiX[dir] * xsiX[dir] - speedOfSoundSq;
		hermSndYY = xsiY[dir] * xsiY[dir] - speedOfSoundSq;
		hermSndZZ = xsiZ[dir] * xsiZ[dir] - speedOfSoundSq;
		hermSndXY = xsiX[dir] * xsiY[dir];
		hermSndXZ = xsiX[dir] * xsiZ[dir];
		hermSndYZ = xsiY[dir] * xsiZ[dir];

		//!calculate third-order hermite tensor elements
		hermThrdXXY = hermSndXX * xsiY[dir];
		hermThrdXYY = hermSndYY * xsiX[dir];
		hermThrdXXZ = hermSndXX * xsiZ[dir];
		hermThrdXZZ = hermSndZZ * xsiX[dir];
		hermThrdYYZ = hermSndYY * xsiZ[dir];
		hermThrdYZZ = hermSndZZ * xsiY[dir];

		//!orthogonalize third-order hermite polynomials between each other
		hermThrdXXYpYZZ = hermThrdXXY + hermThrdYZZ;
		hermThrdXXZpYYZ = hermThrdXXZ + hermThrdYYZ;
		hermThrdXYYpXZZ = hermThrdXYY + hermThrdXZZ;
		hermThrdYZZmXXY = hermThrdYZZ - hermThrdXXY;
		hermThrdYYZmXXZ = hermThrdYYZ - hermThrdXXZ;
		hermThrdXZZmXYY = hermThrdXZZ - hermThrdXYY;

		nEquilHRR[dir] = weight[dir] * (oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aHybridXX + hermSndYY * aHybridYY + hermSndZZ * aHybridZZ
			+ hermSndXY * aHybridXY * 2. + hermSndXZ * aHybridXZ * 2. + hermSndYZ * aHybridYZ * 2.) //!second-order terms from hermite polynomial approximation

						//+ oneOver6 * reciSpeedOfSoundHc * (hermThrdXXY * aNeqThrdXXY * 3.
						  //							   + hermThrdXYY * aNeqThrdXYY * 3.
						  //							   + hermThrdXXZ * aNeqThrdXXZ * 3.
						  //							   + hermThrdXZZ * aNeqThrdXZZ * 3.
						  //						       + hermThrdYYZ * aNeqThrdYYZ * 3.
						  //							   + hermThrdYZZ * aNeqThrdYZZ * 3.) //!third-order terms from hermite polynomial approximation

						//!to avoid spurious couplings and thus stability issues for the D3Q19-GH formalism rely on fully orthogonal basis
						//!orthogonalize fourth-order hermite polynomials with second-order ones and further orthogonalize them between each other
						//see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. G17 - G19

						//!orthogonal third-order hermite polynomials (see Jacob et al. 2018, Journal of Turbulence)
						+oneOver6 * reciSpeedOfSoundHc * (hermThrdXXYpYZZ * aNeqThrdXXYpYZZ * 3.
							+ hermThrdXXZpYYZ * aNeqThrdXXZpYYZ * 3.
							+ hermThrdXYYpXZZ * aNeqThrdXYYpXZZ * 3.
							+ hermThrdYZZmXXY * aNeqThrdYZZmXXY
							+ hermThrdYYZmXXZ * aNeqThrdYYZmXXZ
							+ hermThrdXZZmXYY * aNeqThrdXZZmXYY));
	}
}



void Node::calcHRRoffEquilibrium4_GHbasis(bool& alternating, double* nEquilHRR, double* equil, Control& bc, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared
	double  speedOfSoundQc = pow(bc.getSpeedOfSound(), 4.); //!speed of sound, quartic
	double  speedOfSoundHc = pow(bc.getSpeedOfSound(), 6.); //!speed of sound, hexic
	double  speedOfSoundOc = pow(bc.getSpeedOfSound(), 8.); //!speed of sound, octic
	double  reciSpeedOfSoundSq = 1. / speedOfSoundSq; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundQc = 1. / speedOfSoundQc; //!reciprocal value of speed of sound, quartic
	double	reciSpeedOfSoundHc = 1. / speedOfSoundHc; //!reciprocal value of speed of sound, hexic
	double  reciSpeedOfSoundOc = 1. / speedOfSoundOc; //!reciprocal value of speed of sound, octic
	double  oneOver2 = 1. / 2.;
	double  oneOver4 = 1. / 4.;
	double  oneOver6 = 1. / 6.;
	double  twoOver7 = 2. / 7.;
	double  twoOver5 = 2. / 5.;
	double  eightOver7 = 8. / 7.;
	double  fiftysixOver45 = 56. / 45.;
	double  fortyOver27 = 40. / 27.;
	double  hybridPar = bc.getHybridParameter();
	double	timeStep = bc.getTimeStep() / pow(2., level_);

	//!calculate density at node
	double dens = this->calcDensity(bc, alternating);

	//!calculate dynamic viscosity and inverse pre-factor
	double kinVisc = bc.getKinematicViscosity();
	double dynVisc = kinVisc * dens;
	double invPreFactor = -((timeStep * speedOfSoundSq) / (colFreq_ * kinVisc));

	//!ATTENTION: FIRST ORDER TERMS OF NON-EQUILIBRIUM PROJECTION TO HERMITE-BASIS ARE NEEDED WITH GUO FORCING! (SEE: 10.1016/J.compfluid.2019.06.030)
	//!elements of first-order hermite tensor
	double hermFrstX[19] = { 0. }, hermFrstY[19] = { 0. }, hermFrstZ[19] = { 0. };

	//!elements of second-order hermite tensor
	double hermSndXX, hermSndYY, hermSndZZ, hermSndXY, hermSndXZ, hermSndYZ;

	//!elements of third-order hermite tensor
	double hermThrdXXY, hermThrdXYY, hermThrdXXZ, hermThrdXZZ, hermThrdYYZ, hermThrdYZZ, hermThrdXXYpYZZ, hermThrdXXZpYYZ, hermThrdXYYpXZZ, hermThrdYZZmXXY, hermThrdYYZmXXZ, hermThrdXZZmXYY;

	//!elements of fourth-order hermite tensor (complete projection not possible due to insufficient quadrature using D3Q19)
	double hermFrthXXYY, hermFrthXXZZ, hermFrthYYZZ, hermFrthXXYY4o2, hermFrthXXZZ4o2, hermFrthYYZZ4o2, hermFrthXXYYFO, hermFrthXXZZFO, hermFrthYYZZFO;

	//!hybrid coefficients including finite-difference approximated strain-rate tensor
	double aHybridXX = hybridPar * Pi1Tensor_[0][0] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][0]);
	double aHybridYY = hybridPar * Pi1Tensor_[1][1] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[1][1]);
	double aHybridZZ = hybridPar * Pi1Tensor_[2][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[2][2]);
	double aHybridXY = hybridPar * Pi1Tensor_[0][1] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][1]);
	double aHybridXZ = hybridPar * Pi1Tensor_[0][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[0][2]);
	double aHybridYZ = hybridPar * Pi1Tensor_[1][2] + (1. - hybridPar) * invPreFactor * (2. * dynVisc * strainRateTensorFD_[1][2]);

	//!coefficients of first-order non-equilibrium projection to third-order hermite basis
	double aNeqThrdXXY = 2. * velo[0] * aHybridXY + velo[1] * aHybridXX;
	double aNeqThrdXYY = 2. * velo[1] * aHybridXY + velo[0] * aHybridYY;
	double aNeqThrdXXZ = 2. * velo[0] * aHybridXZ + velo[2] * aHybridXX;
	double aNeqThrdXZZ = 2. * velo[2] * aHybridXZ + velo[0] * aHybridZZ;
	double aNeqThrdYYZ = 2. * velo[1] * aHybridYZ + velo[2] * aHybridYY;
	double aNeqThrdYZZ = 2. * velo[2] * aHybridYZ + velo[1] * aHybridZZ;

	//!orthogonalize third-order hermite polynomials between each other
	double aNeqThrdXXYpYZZ = aNeqThrdXXY + aNeqThrdYZZ;
	double aNeqThrdXXZpYYZ = aNeqThrdXXZ + aNeqThrdYYZ;
	double aNeqThrdXYYpXZZ = aNeqThrdXYY + aNeqThrdXZZ;
	double aNeqThrdYZZmXXY = aNeqThrdYZZ - aNeqThrdXXY;
	double aNeqThrdYYZmXXZ = aNeqThrdYYZ - aNeqThrdXXZ;
	double aNeqThrdXZZmXYY = aNeqThrdXZZ - aNeqThrdXYY;

	//!coefficients of first-order non-equilibrium projection to fourth-order hermite basis via their recursive properties
	double aNeqFrthXXYY = velo[1] * velo[1] * aHybridXX + velo[0] * velo[0] * aHybridYY + 4. * velo[0] * velo[1] * aHybridXY;
	double aNeqFrthXXZZ = velo[2] * velo[2] * aHybridXX + velo[0] * velo[0] * aHybridZZ + 4. * velo[0] * velo[2] * aHybridXZ;
	double aNeqFrthYYZZ = velo[1] * velo[1] * aHybridZZ + velo[2] * velo[2] * aHybridYY + 4. * velo[1] * velo[2] * aHybridYZ;

	//!order 4 orthogonal with respect to order 2
	double aNeqFrthXXYY4o2 = aNeqFrthXXYY + oneOver6 * aHybridZZ;
	double aNeqFrthXXZZ4o2 = aNeqFrthXXZZ + oneOver6 * aHybridYY;
	double aNeqFrthYYZZ4o2 = aNeqFrthYYZZ + oneOver6 * aHybridXX;

    //!order 4 fully orthogonal
	double aNeqFrthXXYYFO = aNeqFrthXXYY4o2;
	double aNeqFrthXXZZFO = aNeqFrthXXZZ4o2 + twoOver7 * aNeqFrthXXYYFO;
	double aNeqFrthYYZZFO = aNeqFrthYYZZ4o2 + twoOver7 * aNeqFrthXXYYFO + twoOver5 * aNeqFrthXXZZFO;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate second-order hermite tensor elements
		hermSndXX = xsiX[dir] * xsiX[dir] - speedOfSoundSq;
		hermSndYY = xsiY[dir] * xsiY[dir] - speedOfSoundSq;
		hermSndZZ = xsiZ[dir] * xsiZ[dir] - speedOfSoundSq;
		hermSndXY = xsiX[dir] * xsiY[dir];
		hermSndXZ = xsiX[dir] * xsiZ[dir];
		hermSndYZ = xsiY[dir] * xsiZ[dir];

		//!calculate third-order hermite tensor elements
		hermThrdXXY = hermSndXX * xsiY[dir];
		hermThrdXYY = hermSndYY * xsiX[dir];
		hermThrdXXZ = hermSndXX * xsiZ[dir];
		hermThrdXZZ = hermSndZZ * xsiX[dir];
		hermThrdYYZ = hermSndYY * xsiZ[dir];
		hermThrdYZZ = hermSndZZ * xsiY[dir];

		//!orthogonalize third-order hermite polynomials between each other
		hermThrdXXYpYZZ = hermThrdXXY + hermThrdYZZ;
		hermThrdXXZpYYZ = hermThrdXXZ + hermThrdYYZ;
		hermThrdXYYpXZZ = hermThrdXYY + hermThrdXZZ;
		hermThrdYZZmXXY = hermThrdYZZ - hermThrdXXY;
		hermThrdYYZmXXZ = hermThrdYYZ - hermThrdXXZ;
		hermThrdXZZmXYY = hermThrdXZZ - hermThrdXYY;

		//!calculate fourth-order hermite tensor elements
		hermFrthXXYY = hermSndXX * hermSndYY;
		hermFrthXXZZ = hermSndXX * hermSndZZ;
		hermFrthYYZZ = hermSndYY * hermSndZZ;

		//!order 4 orthogonal with respect to order 2
		hermFrthXXYY4o2 = hermFrthXXYY + oneOver6 * hermSndZZ;
		hermFrthXXZZ4o2 = hermFrthXXZZ + oneOver6 * hermSndYY;
		hermFrthYYZZ4o2 = hermFrthYYZZ + oneOver6 * hermSndXX;

		//!order 4 fully orthogonal
		hermFrthXXYYFO = hermFrthXXYY4o2;
		hermFrthXXZZFO = hermFrthXXZZ4o2 + twoOver7 * hermFrthXXYYFO;
		hermFrthYYZZFO = hermFrthYYZZ4o2 + twoOver7 * hermFrthXXYYFO + twoOver5 * hermFrthXXZZFO;

		nEquilHRR[dir] = weight[dir] * (oneOver2 * reciSpeedOfSoundQc * (hermSndXX * aHybridXX + hermSndYY * aHybridYY + hermSndZZ * aHybridZZ
			+ hermSndXY * aHybridXY * 2. + hermSndXZ * aHybridXZ * 2. + hermSndYZ * aHybridYZ * 2.) //!second-order terms from hermite polynomial approximation

						//+ oneOver6 * reciSpeedOfSoundHc * (hermThrdXXY * aNeqThrdXXY * 3.
						  //							   + hermThrdXYY * aNeqThrdXYY * 3.
						  //							   + hermThrdXXZ * aNeqThrdXXZ * 3.
						  //							   + hermThrdXZZ * aNeqThrdXZZ * 3.
						  //						       + hermThrdYYZ * aNeqThrdYYZ * 3.
						  //							   + hermThrdYZZ * aNeqThrdYZZ * 3.) //!third-order terms from hermite polynomial approximation
						//+ oneOver4 * reciSpeedOfSoundOc * (hermFrthXXYY * aNeqFrthXXYY + hermFrthXXZZ * aNeqFrthXXZZ + hermFrthYYZZ * aNeqFrthYYZZ)); //!fourth-order terms from hermite polynomial approximation


						//!to avoid spurious couplings and thus stability issues for the D3Q19-GH formalism rely on fully orthogonal basis
						//!orthogonalize fourth-order hermite polynomials with second-order ones and further orthogonalize them between each other
						//see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. G17 - G19

						//!orthogonal third-order hermite polynomials (see Jacob et al. 2018, Journal of Turbulence)
						+oneOver6 * reciSpeedOfSoundHc * (hermThrdXXYpYZZ * aNeqThrdXXYpYZZ * 3.
							+ hermThrdXXZpYYZ * aNeqThrdXXZpYYZ * 3.
							+ hermThrdXYYpXZZ * aNeqThrdXYYpXZZ * 3.
							+ hermThrdYZZmXXY * aNeqThrdYZZmXXY
							+ hermThrdYYZmXXZ * aNeqThrdYYZmXXZ
							+ hermThrdXZZmXYY * aNeqThrdXZZmXYY)

		//!fully orthogonal fourth-order hermite polynomials
		+ oneOver4 * reciSpeedOfSoundOc * (eightOver7 * hermFrthXXYYFO * aNeqFrthXXYYFO + fiftysixOver45 * hermFrthXXZZFO * aNeqFrthXXZZFO + fortyOver27 * hermFrthYYZZFO * aNeqFrthYYZZFO));
	}
}


void Node::calcStrainRateTensorFDandCubicMachCorrection(Control& bc)
{
	double  spacing = bc.getSpacing() / pow(2., level_);
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  cs2 = pow(bc.getSpeedOfSound(), 2.); //!reciprocal value of speed of sound, squared
	double  rec_cs4 = 1. / (cs2 * cs2); //!speed of sound, quartic

	//!elements of second-order hermite tensor
	double hermSndXX[19], hermSndYY[19], hermSndZZ[19], hermSndXY[19], hermSndXZ[19], hermSndYZ[19];

	//!cubic Mach deviation terms at neighbors
	double PsiXXX[6] = { 0. }, PsiXXY[6] = { 0. }, PsiXXZ[6] = { 0. };
	double PsiXYY[6] = { 0. }, PsiYYY[6] = { 0. }, PsiYYZ[6] = { 0. };
	double PsiXZZ[6] = { 0. }, PsiYZZ[6] = { 0. }, PsiZZZ[6] = { 0. };
	double PsiXYZ[6] = { 0. };

	//!spatial derivatives of cubic Mach deviation terms
	double PsiXXX_dx = 0., PsiXXY_dy = 0., PsiXXZ_dz = 0.;
	double PsiXYY_dx = 0., PsiYYY_dy = 0., PsiYYZ_dz = 0.;
	double PsiXZZ_dx = 0., PsiYZZ_dy = 0., PsiZZZ_dz = 0.;
	double PsiXXY_dx = 0., PsiXYY_dy = 0., PsiXYZ_dz = 0.;
	double PsiXXZ_dx = 0., PsiXYZ_dy = 0., PsiXZZ_dz = 0.;
	double PsiXYZ_dx = 0., PsiYYZ_dy = 0., PsiYZZ_dz = 0.;

	//!hermite tensor elements
	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate main diagonal second-order hermite tensor elements
		hermSndXX[dir] = pow(xsiX[dir], 2.) - cs2;
		hermSndYY[dir] = pow(xsiY[dir], 2.) - cs2;
		hermSndZZ[dir] = pow(xsiZ[dir], 2.) - cs2;

		//!calculate secondary diagonal second-order hermite tensor elements
		hermSndXY[dir] = xsiX[dir] * xsiY[dir];
		hermSndXZ[dir] = xsiX[dir] * xsiZ[dir];
		hermSndYZ[dir] = xsiY[dir] * xsiZ[dir];
	}

	//!velocity vectors from direct neighbors in cartesian directions
	double velocityNeighbors[6][3] = { 0. };
	double distanceToNeighbor[6] = { 0. };

	if (this->getTag() == "gridcoupling") //!if current node is instance of gridcoupling
	{
		return; //!since gridcoupling nodes receive their populations and strain-rate values from partner --> return
	}
	else
	{
		//!get neighbor velocities
		for (int dir = 0; dir <= 5; dir++)
		{
			if (this->getNeighbor(dir)->getTag() == "wall")
			{
				//!central difference scheme with assumed zero-velocity at wall!
				//!for wall neighbors set velocity to zero
				velocityNeighbors[dir][0] = 0.;
				velocityNeighbors[dir][1] = 0.;
				velocityNeighbors[dir][2] = 0.;;

				//!calculate cubic Mach devation terms at neighbors
				PsiXXX[dir] = 0.;
				PsiXXY[dir] = 0.;
				PsiXXZ[dir] = 0.;
				PsiXYY[dir] = 0.;
				PsiYYY[dir] = 0.;
				PsiYYZ[dir] = 0.;
				PsiXZZ[dir] = 0.;
				PsiYZZ[dir] = 0.;
				PsiZZZ[dir] = 0.;
				PsiXYZ[dir] = 0.;

				distanceToNeighbor[dir] = 0.5 * spacing;
			}
			else
			{
				double density = this->getNeighbor(dir)->calcDensityPreCol(bc);
				this->getNeighbor(dir)->calcVelocityPreCol(bc, density, velocityNeighbors[dir]);

				//!calculate cubic Mach devation terms at neighbors
				PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
				PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
				PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
				PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
				PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
				PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
				PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
				PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
				PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
				PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];

				distanceToNeighbor[dir] = spacing;
			}
		}
	}

	//!calculate spatial dervatives of cubic Mach corrections terms
	PsiXXX_dx = (PsiXXX[0] - PsiXXX[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYY_dx = (PsiXYY[0] - PsiXYY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXZZ_dx = (PsiXZZ[0] - PsiXZZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXY_dx = (PsiXXY[0] - PsiXXY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXZ_dx = (PsiXXZ[0] - PsiXXZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYZ_dx = (PsiXYZ[0] - PsiXYZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);

	PsiXXY_dy = (PsiXXY[1] - PsiXXY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYY_dy = (PsiYYY[1] - PsiYYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYZZ_dy = (PsiYZZ[1] - PsiYZZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYY_dy = (PsiXYY[1] - PsiXYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYZ_dy = (PsiXYZ[1] - PsiXYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYZ_dy = (PsiYYZ[1] - PsiYYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);

	PsiXXZ_dz = (PsiXXZ[4] - PsiXXZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYYZ_dz = (PsiYYZ[4] - PsiYYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiZZZ_dz = (PsiZZZ[4] - PsiZZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXYZ_dz = (PsiXYZ[4] - PsiXYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXZZ_dz = (PsiXZZ[4] - PsiXZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYZZ_dz = (PsiYZZ[4] - PsiYZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	//!calculate complete cubic Mach correction for all populations
	for (int dir = 0; dir < 19; dir++)
	{
		psi_[dir] = 0.5 * weight[dir] * rec_cs4 * (hermSndXX[dir] * (PsiXXX_dx + PsiXXY_dy + PsiXXZ_dz)
			+ hermSndYY[dir] * (PsiXYY_dx + PsiYYY_dy + PsiYYZ_dz)
			+ hermSndZZ[dir] * (PsiXZZ_dx + PsiYZZ_dy + PsiZZZ_dz)
			+ hermSndXY[dir] * (PsiXXY_dx + PsiXYY_dy + PsiXYZ_dz) * 2.
			+ hermSndXZ[dir] * (PsiXXZ_dx + PsiXYZ_dy + PsiXZZ_dz) * 2.
			+ hermSndYZ[dir] * (PsiXYZ_dx + PsiYYZ_dy + PsiYZZ_dz) * 2.);
	}

	//!calculate strain-rate tensor coefficients via finite difference approximation
	strainRateTensorFD_[0][0] = (velocityNeighbors[0][0] - velocityNeighbors[2][0]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	strainRateTensorFD_[1][1] = (velocityNeighbors[1][1] - velocityNeighbors[3][1]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	strainRateTensorFD_[2][2] = (velocityNeighbors[4][2] - velocityNeighbors[5][2]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	strainRateTensorFD_[0][1] = 0.5 * ((velocityNeighbors[0][1] - velocityNeighbors[2][1]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[1][0] - velocityNeighbors[3][0]) / (distanceToNeighbor[1] + distanceToNeighbor[3]));
	strainRateTensorFD_[0][2] = 0.5 * ((velocityNeighbors[0][2] - velocityNeighbors[2][2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[4][0] - velocityNeighbors[5][0]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
	strainRateTensorFD_[1][2] = 0.5 * ((velocityNeighbors[1][2] - velocityNeighbors[3][2]) / (distanceToNeighbor[1] + distanceToNeighbor[3]) + (velocityNeighbors[4][1] - velocityNeighbors[5][1]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
}


void Node::calcStrainRateTensorFDandCubicMachCorrectionDirectCoupling(bool& alternatingL0, bool& alternatingL1, Control& bc)
{
	double  spacing = bc.getSpacing() / pow(2., level_);
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  cs2 = pow(bc.getSpeedOfSound(), 2.); //!reciprocal value of speed of sound, squared
	double  rec_cs4 = 1. / (cs2 * cs2); //!speed of sound, quartic

	//!elements of second-order hermite tensor
	double hermSndXX[19], hermSndYY[19], hermSndZZ[19], hermSndXY[19], hermSndXZ[19], hermSndYZ[19];

	//!cubic Mach deviation terms at neighbors
	double PsiXXX[6] = { 0. }, PsiXXY[6] = { 0. }, PsiXXZ[6] = { 0. };
	double PsiXYY[6] = { 0. }, PsiYYY[6] = { 0. }, PsiYYZ[6] = { 0. };
	double PsiXZZ[6] = { 0. }, PsiYZZ[6] = { 0. }, PsiZZZ[6] = { 0. };
	double PsiXYZ[6] = { 0. };

	//!spatial derivatives of cubic Mach deviation terms
	double PsiXXX_dx = 0., PsiXXY_dy = 0., PsiXXZ_dz = 0.;
	double PsiXYY_dx = 0., PsiYYY_dy = 0., PsiYYZ_dz = 0.;
	double PsiXZZ_dx = 0., PsiYZZ_dy = 0., PsiZZZ_dz = 0.;
	double PsiXXY_dx = 0., PsiXYY_dy = 0., PsiXYZ_dz = 0.;
	double PsiXXZ_dx = 0., PsiXYZ_dy = 0., PsiXZZ_dz = 0.;
	double PsiXYZ_dx = 0., PsiYYZ_dy = 0., PsiYZZ_dz = 0.;

	//!hermite tensor elements
	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate main diagonal second-order hermite tensor elements
		hermSndXX[dir] = pow(xsiX[dir], 2.) - cs2;
		hermSndYY[dir] = pow(xsiY[dir], 2.) - cs2;
		hermSndZZ[dir] = pow(xsiZ[dir], 2.) - cs2;

		//!calculate secondary diagonal second-order hermite tensor elements
		hermSndXY[dir] = xsiX[dir] * xsiY[dir];
		hermSndXZ[dir] = xsiX[dir] * xsiZ[dir];
		hermSndYZ[dir] = xsiY[dir] * xsiZ[dir];
	}

	//!velocity vectors from direct neighbors in cartesian directions
	double velocityNeighbors[6][3] = { 0. };
	double distanceToNeighbor[6] = { 0. };

	//!get neighbor velocities
	for (int dir = 0; dir <= 5; dir++)
	{
		if ((this->getNeighbor(dir) == NULL) && (level_ == 0))
		{
			//!if coarse neighbor doesn't exist, change level and get temporally interpolated value of fine partner's dir->dir neighbor node
			double density = this->getPartner()->getNeighbor(dir)->getNeighbor(dir)->calcDensity(bc, alternatingL1);
			this->getPartner()->getNeighbor(dir)->getNeighbor(dir)->calcVelocity(bc, alternatingL1, density, velocityNeighbors[dir]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];

			distanceToNeighbor[dir] = spacing;
		}
		else if ((this->getNeighbor(dir) != NULL) && (level_ == 0))
		{
			double density = this->getNeighbor(dir)->calcDensity(bc, alternatingL0);
			this->getNeighbor(dir)->calcVelocity(bc, alternatingL0, density, velocityNeighbors[dir]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];

			distanceToNeighbor[dir] = spacing;
		}
		else if ((this->getNeighbor(dir) == NULL) && (level_ == 1))
		{
			//!this routine is called on fine interface nodes only during asynchronous iteration
			//!if fine neighbor doesn't exist, change level and get temporally interpolated value of coarse partner's dir neighbor node
			//!during synchronous iteration, the strain-rate tensor and correction terms are calculated on the coarse node and transferred to it's fine partner (see below)
			double veloPast[3] = { 0. };
			double densityPast = this->getPartner()->getNeighbor(dir)->calcDensityPreCol(bc);
			this->getPartner()->getNeighbor(dir)->calcVelocityPreCol(bc, densityPast, veloPast);

			double veloFuture[3] = { 0. };
			double densityFuture = this->getPartner()->getNeighbor(dir)->calcDensity(bc, alternatingL0);
			this->getPartner()->getNeighbor(dir)->calcVelocity(bc, alternatingL0, densityFuture, veloFuture);

			double density = 0.5 * (densityPast + densityFuture);
			velocityNeighbors[dir][0] = 0.5 * (veloPast[0] + veloFuture[0]);
			velocityNeighbors[dir][1] = 0.5 * (veloPast[1] + veloFuture[1]);
			velocityNeighbors[dir][2] = 0.5 * (veloPast[2] + veloFuture[2]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];

			//!centered finite differences need to be calculated with coarse stencil always!
			distanceToNeighbor[dir] = 2. * spacing;
		}
		else if ((this->getNeighbor(dir) != NULL) && (level_ == 1))
		{
			double density = this->getNeighbor(dir)->getNeighbor(dir)->calcDensity(bc, alternatingL1);
			this->getNeighbor(dir)->getNeighbor(dir)->calcVelocity(bc, alternatingL1, density, velocityNeighbors[dir]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];

			//!centered finite differences need to be calculated with coarse stencil always!
			distanceToNeighbor[dir] = 2. * spacing;
		}
	}

	//!calculate spatial dervatives of cubic Mach corrections terms
	PsiXXX_dx = (PsiXXX[0] - PsiXXX[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYY_dx = (PsiXYY[0] - PsiXYY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXZZ_dx = (PsiXZZ[0] - PsiXZZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXY_dx = (PsiXXY[0] - PsiXXY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXZ_dx = (PsiXXZ[0] - PsiXXZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYZ_dx = (PsiXYZ[0] - PsiXYZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);

	PsiXXY_dy = (PsiXXY[1] - PsiXXY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYY_dy = (PsiYYY[1] - PsiYYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYZZ_dy = (PsiYZZ[1] - PsiYZZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYY_dy = (PsiXYY[1] - PsiXYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYZ_dy = (PsiXYZ[1] - PsiXYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYZ_dy = (PsiYYZ[1] - PsiYYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);

	PsiXXZ_dz = (PsiXXZ[4] - PsiXXZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYYZ_dz = (PsiYYZ[4] - PsiYYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiZZZ_dz = (PsiZZZ[4] - PsiZZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXYZ_dz = (PsiXYZ[4] - PsiXYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXZZ_dz = (PsiXZZ[4] - PsiXZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYZZ_dz = (PsiYZZ[4] - PsiYZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	//!if function is called on coarse interface node --> synchronous iteration --> transfer strain-rates and correction terms to fine partner
	if (level_ == 0)
	{
		//!calculate complete cubic Mach correction for all populations
		for (int dir = 0; dir < 19; dir++)
		{
			psi_[dir] = 0.5 * weight[dir] * rec_cs4 * (hermSndXX[dir] * (PsiXXX_dx + PsiXXY_dy + PsiXXZ_dz)
													 + hermSndYY[dir] * (PsiXYY_dx + PsiYYY_dy + PsiYYZ_dz)
													 + hermSndZZ[dir] * (PsiXZZ_dx + PsiYZZ_dy + PsiZZZ_dz)
													 + hermSndXY[dir] * (PsiXXY_dx + PsiXYY_dy + PsiXYZ_dz) * 2.
													 + hermSndXZ[dir] * (PsiXXZ_dx + PsiXYZ_dy + PsiXZZ_dz) * 2.
													 + hermSndYZ[dir] * (PsiXYZ_dx + PsiYYZ_dy + PsiYZZ_dz) * 2.);

			this->getPartner()->setCubicMachCorrection(dir, psi_[dir]);
		}

		//!calculate strain-rate tensor coefficients via finite difference approximation
		strainRateTensorFD_[0][0] = (velocityNeighbors[0][0] - velocityNeighbors[2][0]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
		strainRateTensorFD_[1][1] = (velocityNeighbors[1][1] - velocityNeighbors[3][1]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
		strainRateTensorFD_[2][2] = (velocityNeighbors[4][2] - velocityNeighbors[5][2]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

		strainRateTensorFD_[0][1] = 0.5 * ((velocityNeighbors[0][1] - velocityNeighbors[2][1]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[1][0] - velocityNeighbors[3][0]) / (distanceToNeighbor[1] + distanceToNeighbor[3]));
		strainRateTensorFD_[0][2] = 0.5 * ((velocityNeighbors[0][2] - velocityNeighbors[2][2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[4][0] - velocityNeighbors[5][0]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
		strainRateTensorFD_[1][2] = 0.5 * ((velocityNeighbors[1][2] - velocityNeighbors[3][2]) / (distanceToNeighbor[1] + distanceToNeighbor[3]) + (velocityNeighbors[4][1] - velocityNeighbors[5][1]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));

		this->getPartner()->setStrainRateTensorFD(0, 0, strainRateTensorFD_[0][0]);
		this->getPartner()->setStrainRateTensorFD(1, 1, strainRateTensorFD_[1][1]);
		this->getPartner()->setStrainRateTensorFD(2, 2, strainRateTensorFD_[2][2]);

		this->getPartner()->setStrainRateTensorFD(0, 1, strainRateTensorFD_[0][1]);
		this->getPartner()->setStrainRateTensorFD(0, 2, strainRateTensorFD_[0][2]);
		this->getPartner()->setStrainRateTensorFD(1, 2, strainRateTensorFD_[1][2]);
	}
	else
	{
		//!calculate complete cubic Mach correction for all populations
		for (int dir = 0; dir < 19; dir++)
		{
			psi_[dir] = 0.5 * weight[dir] * rec_cs4 * (hermSndXX[dir] * (PsiXXX_dx + PsiXXY_dy + PsiXXZ_dz)
													 + hermSndYY[dir] * (PsiXYY_dx + PsiYYY_dy + PsiYYZ_dz)
													 + hermSndZZ[dir] * (PsiXZZ_dx + PsiYZZ_dy + PsiZZZ_dz)
													 + hermSndXY[dir] * (PsiXXY_dx + PsiXYY_dy + PsiXYZ_dz) * 2.
													 + hermSndXZ[dir] * (PsiXXZ_dx + PsiXYZ_dy + PsiXZZ_dz) * 2.
													 + hermSndYZ[dir] * (PsiXYZ_dx + PsiYYZ_dy + PsiYZZ_dz) * 2.);
		}

		//!calculate strain-rate tensor coefficients via finite difference approximation
		strainRateTensorFD_[0][0] = (velocityNeighbors[0][0] - velocityNeighbors[2][0]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
		strainRateTensorFD_[1][1] = (velocityNeighbors[1][1] - velocityNeighbors[3][1]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
		strainRateTensorFD_[2][2] = (velocityNeighbors[4][2] - velocityNeighbors[5][2]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

		strainRateTensorFD_[0][1] = 0.5 * ((velocityNeighbors[0][1] - velocityNeighbors[2][1]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[1][0] - velocityNeighbors[3][0]) / (distanceToNeighbor[1] + distanceToNeighbor[3]));
		strainRateTensorFD_[0][2] = 0.5 * ((velocityNeighbors[0][2] - velocityNeighbors[2][2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[4][0] - velocityNeighbors[5][0]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
		strainRateTensorFD_[1][2] = 0.5 * ((velocityNeighbors[1][2] - velocityNeighbors[3][2]) / (distanceToNeighbor[1] + distanceToNeighbor[3]) + (velocityNeighbors[4][1] - velocityNeighbors[5][1]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
	}
}


void Node::calcStrainRateTensorFDandCubicMachCorrectionCoarsePartners(bool& alternatingL0, bool& alternatingL1, Control& bc)
{
	double  spacing = bc.getSpacing() / pow(2., level_);
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  cs2 = pow(bc.getSpeedOfSound(), 2.); //!reciprocal value of speed of sound, squared
	double  rec_cs4 = 1. / (cs2 * cs2); //!speed of sound, quartic

	//!elements of second-order hermite tensor
	double hermSndXX[19], hermSndYY[19], hermSndZZ[19], hermSndXY[19], hermSndXZ[19], hermSndYZ[19];

	//!cubic Mach deviation terms at neighbors
	double PsiXXX[6] = { 0. }, PsiXXY[6] = { 0. }, PsiXXZ[6] = { 0. };
	double PsiXYY[6] = { 0. }, PsiYYY[6] = { 0. }, PsiYYZ[6] = { 0. };
	double PsiXZZ[6] = { 0. }, PsiYZZ[6] = { 0. }, PsiZZZ[6] = { 0. };
	double PsiXYZ[6] = { 0. };

	//!spatial derivatives of cubic Mach deviation terms
	double PsiXXX_dx = 0., PsiXXY_dy = 0., PsiXXZ_dz = 0.;
	double PsiXYY_dx = 0., PsiYYY_dy = 0., PsiYYZ_dz = 0.;
	double PsiXZZ_dx = 0., PsiYZZ_dy = 0., PsiZZZ_dz = 0.;
	double PsiXXY_dx = 0., PsiXYY_dy = 0., PsiXYZ_dz = 0.;
	double PsiXXZ_dx = 0., PsiXYZ_dy = 0., PsiXZZ_dz = 0.;
	double PsiXYZ_dx = 0., PsiYYZ_dy = 0., PsiYZZ_dz = 0.;

	//!hermite tensor elements
	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate main diagonal second-order hermite tensor elements
		hermSndXX[dir] = pow(xsiX[dir], 2.) - cs2;
		hermSndYY[dir] = pow(xsiY[dir], 2.) - cs2;
		hermSndZZ[dir] = pow(xsiZ[dir], 2.) - cs2;

		//!calculate secondary diagonal second-order hermite tensor elements
		hermSndXY[dir] = xsiX[dir] * xsiY[dir];
		hermSndXZ[dir] = xsiX[dir] * xsiZ[dir];
		hermSndYZ[dir] = xsiY[dir] * xsiZ[dir];
	}

	//!velocity vectors from direct neighbors in cartesian directions
	double velocityNeighbors[6][3] = { 0. };
	double distanceToNeighbor[6] = { 0. };

	//!get neighbor velocities
	for (int dir = 0; dir <= 5; dir++)
	{
		if (this->getNeighbor(dir)->getTag() == "gridcoupling")
		{
			double density = this->getNeighbor(dir)->getPartner()->calcDensity(bc, alternatingL1);
			this->getNeighbor(dir)->getPartner()->calcVelocity(bc, alternatingL1, density, velocityNeighbors[dir]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];
		}
		else
		{
			double density = this->getNeighbor(dir)->calcDensity(bc, alternatingL0);
			this->getNeighbor(dir)->calcVelocity(bc, alternatingL0, density, velocityNeighbors[dir]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];
		}
		distanceToNeighbor[dir] = spacing;
	}

	//!calculate spatial dervatives of cubic Mach corrections terms
	PsiXXX_dx = (PsiXXX[0] - PsiXXX[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYY_dx = (PsiXYY[0] - PsiXYY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXZZ_dx = (PsiXZZ[0] - PsiXZZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXY_dx = (PsiXXY[0] - PsiXXY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXZ_dx = (PsiXXZ[0] - PsiXXZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYZ_dx = (PsiXYZ[0] - PsiXYZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);

	PsiXXY_dy = (PsiXXY[1] - PsiXXY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYY_dy = (PsiYYY[1] - PsiYYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYZZ_dy = (PsiYZZ[1] - PsiYZZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYY_dy = (PsiXYY[1] - PsiXYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYZ_dy = (PsiXYZ[1] - PsiXYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYZ_dy = (PsiYYZ[1] - PsiYYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);

	PsiXXZ_dz = (PsiXXZ[4] - PsiXXZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYYZ_dz = (PsiYYZ[4] - PsiYYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiZZZ_dz = (PsiZZZ[4] - PsiZZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXYZ_dz = (PsiXYZ[4] - PsiXYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXZZ_dz = (PsiXZZ[4] - PsiXZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYZZ_dz = (PsiYZZ[4] - PsiYZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	//!calculate complete cubic Mach correction for all populations
	for (int dir = 0; dir < 19; dir++)
	{
		psiFull_[dir] = 0.5 * weight[dir] * rec_cs4 * (hermSndXX[dir] * (PsiXXX_dx + PsiXXY_dy + PsiXXZ_dz)
			+ hermSndYY[dir] * (PsiXYY_dx + PsiYYY_dy + PsiYYZ_dz)
			+ hermSndZZ[dir] * (PsiXZZ_dx + PsiYZZ_dy + PsiZZZ_dz)
			+ hermSndXY[dir] * (PsiXXY_dx + PsiXYY_dy + PsiXYZ_dz) * 2.
			+ hermSndXZ[dir] * (PsiXXZ_dx + PsiXYZ_dy + PsiXZZ_dz) * 2.
			+ hermSndYZ[dir] * (PsiXYZ_dx + PsiYYZ_dy + PsiYZZ_dz) * 2.);
	}

	//!calculate strain-rate tensor coefficients via finite difference approximation
	strainRateTensorFDfull_[0][0] = (velocityNeighbors[0][0] - velocityNeighbors[2][0]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	strainRateTensorFDfull_[1][1] = (velocityNeighbors[1][1] - velocityNeighbors[3][1]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	strainRateTensorFDfull_[2][2] = (velocityNeighbors[4][2] - velocityNeighbors[5][2]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	strainRateTensorFDfull_[0][1] = 0.5 * ((velocityNeighbors[0][1] - velocityNeighbors[2][1]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[1][0] - velocityNeighbors[3][0]) / (distanceToNeighbor[1] + distanceToNeighbor[3]));
	strainRateTensorFDfull_[0][2] = 0.5 * ((velocityNeighbors[0][2] - velocityNeighbors[2][2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[4][0] - velocityNeighbors[5][0]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
	strainRateTensorFDfull_[1][2] = 0.5 * ((velocityNeighbors[1][2] - velocityNeighbors[3][2]) / (distanceToNeighbor[1] + distanceToNeighbor[3]) + (velocityNeighbors[4][1] - velocityNeighbors[5][1]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
}



void Node::calcStrainRateTensorFDandCubicMachCorrectionHalfCoarsePartners(bool& alternatingL0, bool& alternatingL1, Control& bc) //!alternative c2f strain-rate-tensor calculation without fictitious streaming
{
	double  spacing = bc.getSpacing() / pow(2., level_);
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  cs2 = pow(bc.getSpeedOfSound(), 2.); //!reciprocal value of speed of sound, squared
	double  rec_cs4 = 1. / (cs2 * cs2); //!speed of sound, quartic

	//!elements of second-order hermite tensor
	double hermSndXX[19], hermSndYY[19], hermSndZZ[19], hermSndXY[19], hermSndXZ[19], hermSndYZ[19];

	//!cubic Mach deviation terms at neighbors
	double PsiXXX[6] = { 0. }, PsiXXY[6] = { 0. }, PsiXXZ[6] = { 0. };
	double PsiXYY[6] = { 0. }, PsiYYY[6] = { 0. }, PsiYYZ[6] = { 0. };
	double PsiXZZ[6] = { 0. }, PsiYZZ[6] = { 0. }, PsiZZZ[6] = { 0. };
	double PsiXYZ[6] = { 0. };

	//!spatial derivatives of cubic Mach deviation terms
	double PsiXXX_dx = 0., PsiXXY_dy = 0., PsiXXZ_dz = 0.;
	double PsiXYY_dx = 0., PsiYYY_dy = 0., PsiYYZ_dz = 0.;
	double PsiXZZ_dx = 0., PsiYZZ_dy = 0., PsiZZZ_dz = 0.;
	double PsiXXY_dx = 0., PsiXYY_dy = 0., PsiXYZ_dz = 0.;
	double PsiXXZ_dx = 0., PsiXYZ_dy = 0., PsiXZZ_dz = 0.;
	double PsiXYZ_dx = 0., PsiYYZ_dy = 0., PsiYZZ_dz = 0.;

	//!hermite tensor elements
	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate main diagonal second-order hermite tensor elements
		hermSndXX[dir] = pow(xsiX[dir], 2.) - cs2;
		hermSndYY[dir] = pow(xsiY[dir], 2.) - cs2;
		hermSndZZ[dir] = pow(xsiZ[dir], 2.) - cs2;

		//!calculate secondary diagonal second-order hermite tensor elements
		hermSndXY[dir] = xsiX[dir] * xsiY[dir];
		hermSndXZ[dir] = xsiX[dir] * xsiZ[dir];
		hermSndYZ[dir] = xsiY[dir] * xsiZ[dir];
	}

	//!velocity vectors from direct neighbors in cartesian directions
	double velocityNeighbors[6][3] = { 0. };
	double distanceToNeighbor[6] = { 0. };

	//!get neighbor velocities
	for (int dir = 0; dir <= 5; dir++)
	{
		if (this->getNeighbor(dir)->getTag() == "gridcoupling") //!if coarse neighbor is interface node --> get t + 1/2 velocity from fine partner
		{
			double density = this->getNeighbor(dir)->getPartner()->calcDensity(bc, alternatingL1);
			this->getNeighbor(dir)->getPartner()->calcVelocity(bc, alternatingL1, density, velocityNeighbors[dir]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];
		}
		else //!if coarse neighbor is not on interface --> perform linear time interpolation of velocity
		{
			double veloPast[3] = { 0. };
			double densityPast = this->getNeighbor(dir)->calcDensityPreCol(bc);
			this->getNeighbor(dir)->calcVelocityPreCol(bc, densityPast, veloPast);

			double veloFuture[3] = { 0. };
			double densityFuture = this->getNeighbor(dir)->calcDensity(bc, alternatingL0);
			this->getNeighbor(dir)->calcVelocity(bc, alternatingL0, densityFuture, veloFuture);

			double density = 0.5 * (densityPast + densityFuture);
			velocityNeighbors[dir][0] = 0.5 * (veloPast[0] + veloFuture[0]);
			velocityNeighbors[dir][1] = 0.5 * (veloPast[1] + veloFuture[1]);
			velocityNeighbors[dir][2] = 0.5 * (veloPast[2] + veloFuture[2]);

			//!calculate cubic Mach devation terms at neighbors
			PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
			PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
			PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
			PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];
		}
		distanceToNeighbor[dir] = spacing;
	}

	//!calculate spatial dervatives of cubic Mach corrections terms
	PsiXXX_dx = (PsiXXX[0] - PsiXXX[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYY_dx = (PsiXYY[0] - PsiXYY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXZZ_dx = (PsiXZZ[0] - PsiXZZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXY_dx = (PsiXXY[0] - PsiXXY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXZ_dx = (PsiXXZ[0] - PsiXXZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYZ_dx = (PsiXYZ[0] - PsiXYZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);

	PsiXXY_dy = (PsiXXY[1] - PsiXXY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYY_dy = (PsiYYY[1] - PsiYYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYZZ_dy = (PsiYZZ[1] - PsiYZZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYY_dy = (PsiXYY[1] - PsiXYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYZ_dy = (PsiXYZ[1] - PsiXYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYZ_dy = (PsiYYZ[1] - PsiYYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);

	PsiXXZ_dz = (PsiXXZ[4] - PsiXXZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYYZ_dz = (PsiYYZ[4] - PsiYYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiZZZ_dz = (PsiZZZ[4] - PsiZZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXYZ_dz = (PsiXYZ[4] - PsiXYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXZZ_dz = (PsiXZZ[4] - PsiXZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYZZ_dz = (PsiYZZ[4] - PsiYZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	//!calculate complete cubic Mach correction for all populations
	for (int dir = 0; dir < 19; dir++)
	{
		psiHalf_[dir] = 0.5 * weight[dir] * rec_cs4 * (hermSndXX[dir] * (PsiXXX_dx + PsiXXY_dy + PsiXXZ_dz)
			+ hermSndYY[dir] * (PsiXYY_dx + PsiYYY_dy + PsiYYZ_dz)
			+ hermSndZZ[dir] * (PsiXZZ_dx + PsiYZZ_dy + PsiZZZ_dz)
			+ hermSndXY[dir] * (PsiXXY_dx + PsiXYY_dy + PsiXYZ_dz) * 2.
			+ hermSndXZ[dir] * (PsiXXZ_dx + PsiXYZ_dy + PsiXZZ_dz) * 2.
			+ hermSndYZ[dir] * (PsiXYZ_dx + PsiYYZ_dy + PsiYZZ_dz) * 2.);
	}

	//!calculate strain-rate tensor coefficients via finite difference approximation
	strainRateTensorFDhalf_[0][0] = (velocityNeighbors[0][0] - velocityNeighbors[2][0]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	strainRateTensorFDhalf_[1][1] = (velocityNeighbors[1][1] - velocityNeighbors[3][1]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	strainRateTensorFDhalf_[2][2] = (velocityNeighbors[4][2] - velocityNeighbors[5][2]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	strainRateTensorFDhalf_[0][1] = 0.5 * ((velocityNeighbors[0][1] - velocityNeighbors[2][1]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[1][0] - velocityNeighbors[3][0]) / (distanceToNeighbor[1] + distanceToNeighbor[3]));
	strainRateTensorFDhalf_[0][2] = 0.5 * ((velocityNeighbors[0][2] - velocityNeighbors[2][2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[4][0] - velocityNeighbors[5][0]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
	strainRateTensorFDhalf_[1][2] = 0.5 * ((velocityNeighbors[1][2] - velocityNeighbors[3][2]) / (distanceToNeighbor[1] + distanceToNeighbor[3]) + (velocityNeighbors[4][1] - velocityNeighbors[5][1]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
}


void Node::calcStrainRateTensorFDandCubicMachCorrectionFinePartners(bool& alternating, Control& bc)
{
	double  spacing = bc.getSpacing() / pow(2., level_);
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  cs2 = pow(bc.getSpeedOfSound(), 2.); //!reciprocal value of speed of sound, squared
	double  rec_cs4 = 1. / (cs2 * cs2); //!speed of sound, quartic

	//!elements of second-order hermite tensor
	double hermSndXX[19], hermSndYY[19], hermSndZZ[19], hermSndXY[19], hermSndXZ[19], hermSndYZ[19];

	//!cubic Mach deviation terms at neighbors
	double PsiXXX[6] = { 0. }, PsiXXY[6] = { 0. }, PsiXXZ[6] = { 0. };
	double PsiXYY[6] = { 0. }, PsiYYY[6] = { 0. }, PsiYYZ[6] = { 0. };
	double PsiXZZ[6] = { 0. }, PsiYZZ[6] = { 0. }, PsiZZZ[6] = { 0. };
	double PsiXYZ[6] = { 0. };

	//!spatial derivatives of cubic Mach deviation terms
	double PsiXXX_dx = 0., PsiXXY_dy = 0., PsiXXZ_dz = 0.;
	double PsiXYY_dx = 0., PsiYYY_dy = 0., PsiYYZ_dz = 0.;
	double PsiXZZ_dx = 0., PsiYZZ_dy = 0., PsiZZZ_dz = 0.;
	double PsiXXY_dx = 0., PsiXYY_dy = 0., PsiXYZ_dz = 0.;
	double PsiXXZ_dx = 0., PsiXYZ_dy = 0., PsiXZZ_dz = 0.;
	double PsiXYZ_dx = 0., PsiYYZ_dy = 0., PsiYZZ_dz = 0.;

	//!hermite tensor elements
	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate main diagonal second-order hermite tensor elements
		hermSndXX[dir] = pow(xsiX[dir], 2.) - cs2;
		hermSndYY[dir] = pow(xsiY[dir], 2.) - cs2;
		hermSndZZ[dir] = pow(xsiZ[dir], 2.) - cs2;

		//!calculate secondary diagonal second-order hermite tensor elements
		hermSndXY[dir] = xsiX[dir] * xsiY[dir];
		hermSndXZ[dir] = xsiX[dir] * xsiZ[dir];
		hermSndYZ[dir] = xsiY[dir] * xsiZ[dir];
	}

	//!velocity vectors from direct neighbors in cartesian directions
	double velocityNeighbors[6][3] = { 0. };
	double distanceToNeighbor[6] = { 0. };

	//!get neighbor velocities
	for (int dir = 0; dir <= 5; dir++)
	{
		double density = this->getNeighbor(dir)->calcDensity(bc, alternating);
		this->getNeighbor(dir)->calcVelocity(bc, alternating, density, velocityNeighbors[dir]);

		//!calculate cubic Mach devation terms at neighbors
		PsiXXX[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
		PsiXXY[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
		PsiXXZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
		PsiXYY[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
		PsiYYY[dir] = density * velocityNeighbors[dir][1] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
		PsiYYZ[dir] = 0.5 * density * velocityNeighbors[dir][2] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
		PsiXZZ[dir] = 0.5 * density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][1];
		PsiYZZ[dir] = 0.5 * density * velocityNeighbors[dir][1] * velocityNeighbors[dir][0] * velocityNeighbors[dir][0];
		PsiZZZ[dir] = density * velocityNeighbors[dir][2] * velocityNeighbors[dir][2] * velocityNeighbors[dir][2];
		PsiXYZ[dir] = density * velocityNeighbors[dir][0] * velocityNeighbors[dir][1] * velocityNeighbors[dir][2];

		distanceToNeighbor[dir] = spacing;
	}

	//!calculate spatial dervatives of cubic Mach corrections terms
	PsiXXX_dx = (PsiXXX[0] - PsiXXX[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYY_dx = (PsiXYY[0] - PsiXYY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXZZ_dx = (PsiXZZ[0] - PsiXZZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXY_dx = (PsiXXY[0] - PsiXXY[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXXZ_dx = (PsiXXZ[0] - PsiXXZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	PsiXYZ_dx = (PsiXYZ[0] - PsiXYZ[2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);

	PsiXXY_dy = (PsiXXY[1] - PsiXXY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYY_dy = (PsiYYY[1] - PsiYYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYZZ_dy = (PsiYZZ[1] - PsiYZZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYY_dy = (PsiXYY[1] - PsiXYY[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiXYZ_dy = (PsiXYZ[1] - PsiXYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	PsiYYZ_dy = (PsiYYZ[1] - PsiYYZ[3]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);

	PsiXXZ_dz = (PsiXXZ[4] - PsiXXZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYYZ_dz = (PsiYYZ[4] - PsiYYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiZZZ_dz = (PsiZZZ[4] - PsiZZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXYZ_dz = (PsiXYZ[4] - PsiXYZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiXZZ_dz = (PsiXZZ[4] - PsiXZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);
	PsiYZZ_dz = (PsiYZZ[4] - PsiYZZ[5]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	//!calculate complete cubic Mach correction for all populations
	for (int dir = 0; dir < 19; dir++)
	{
		psi_[dir] = 0.5 * weight[dir] * rec_cs4 * (hermSndXX[dir] * (PsiXXX_dx + PsiXXY_dy + PsiXXZ_dz)
			+ hermSndYY[dir] * (PsiXYY_dx + PsiYYY_dy + PsiYYZ_dz)
			+ hermSndZZ[dir] * (PsiXZZ_dx + PsiYZZ_dy + PsiZZZ_dz)
			+ hermSndXY[dir] * (PsiXXY_dx + PsiXYY_dy + PsiXYZ_dz) * 2.
			+ hermSndXZ[dir] * (PsiXXZ_dx + PsiXYZ_dy + PsiXZZ_dz) * 2.
			+ hermSndYZ[dir] * (PsiXYZ_dx + PsiYYZ_dy + PsiYZZ_dz) * 2.);
	}

	//!calculate strain-rate tensor coefficients via finite difference approximation
	strainRateTensorFD_[0][0] = (velocityNeighbors[0][0] - velocityNeighbors[2][0]) / (distanceToNeighbor[0] + distanceToNeighbor[2]);
	strainRateTensorFD_[1][1] = (velocityNeighbors[1][1] - velocityNeighbors[3][1]) / (distanceToNeighbor[1] + distanceToNeighbor[3]);
	strainRateTensorFD_[2][2] = (velocityNeighbors[4][2] - velocityNeighbors[5][2]) / (distanceToNeighbor[4] + distanceToNeighbor[5]);

	strainRateTensorFD_[0][1] = 0.5 * ((velocityNeighbors[0][1] - velocityNeighbors[2][1]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[1][0] - velocityNeighbors[3][0]) / (distanceToNeighbor[1] + distanceToNeighbor[3]));
	strainRateTensorFD_[0][2] = 0.5 * ((velocityNeighbors[0][2] - velocityNeighbors[2][2]) / (distanceToNeighbor[0] + distanceToNeighbor[2]) + (velocityNeighbors[4][0] - velocityNeighbors[5][0]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
	strainRateTensorFD_[1][2] = 0.5 * ((velocityNeighbors[1][2] - velocityNeighbors[3][2]) / (distanceToNeighbor[1] + distanceToNeighbor[3]) + (velocityNeighbors[4][1] - velocityNeighbors[5][1]) / (distanceToNeighbor[4] + distanceToNeighbor[5]));
}


void Node::calcPi1Tensor(bool& alternating, double* equil, Control& bc, double* velo)
{
	double* weight = bc.getWeightFactors(); //!lattice weightings
	double* xsiX = bc.getXsiX();
	double* xsiY = bc.getXsiY();
	double* xsiZ = bc.getXsiZ();
	double  speedOfSoundSq = pow(bc.getSpeedOfSound(), 2.); //!speed of sound, squared

	//!elements of second-order hermite tensor
	double hermSndXX[19], hermSndYY[19], hermSndZZ[19], hermSndXY[19], hermSndXZ[19], hermSndYZ[19];

	//!coefficients of first-order off-equilibrium projection to second-order hermite basis
	double aNeqSndXX = 0.;
	double aNeqSndYY = 0.;
	double aNeqSndZZ = 0.;
	double aNeqSndXY = 0.;
	double aNeqSndXZ = 0.;
	double aNeqSndYZ = 0.;

	for (int dir = 0; dir < 19; dir++)
	{
		//!calculate main diagonal second-order hermite tensor elements
		hermSndXX[dir] = pow(xsiX[dir], 2.) - speedOfSoundSq;
		hermSndYY[dir] = pow(xsiY[dir], 2.) - speedOfSoundSq;
		hermSndZZ[dir] = pow(xsiZ[dir], 2.) - speedOfSoundSq;

		//!calculate secondary diagonal second-order hermite tensor elements
		hermSndXY[dir] = xsiX[dir] * xsiY[dir];
		hermSndXZ[dir] = xsiX[dir] * xsiZ[dir];
		hermSndYZ[dir] = xsiY[dir] * xsiZ[dir];

		if (bc.getCubicMachCorrection() == "yes")
		{
			//!calculate main second-order coefficients
			aNeqSndXX += hermSndXX[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil) + 0.5 * timeStep_ * psi_[dir]);
			aNeqSndYY += hermSndYY[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil) + 0.5 * timeStep_ * psi_[dir]);
			aNeqSndZZ += hermSndZZ[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil) + 0.5 * timeStep_ * psi_[dir]);

			//!calculate secondary second-order coefficients
			aNeqSndXY += hermSndXY[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil) + 0.5 * timeStep_ * psi_[dir]);
			aNeqSndXZ += hermSndXZ[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil) + 0.5 * timeStep_ * psi_[dir]);
			aNeqSndYZ += hermSndYZ[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil) + 0.5 * timeStep_ * psi_[dir]);
		}
		else
		{
			//!calculate main second-order coefficients
			aNeqSndXX += hermSndXX[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil));
			aNeqSndYY += hermSndYY[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil));
			aNeqSndZZ += hermSndZZ[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil));

			//!calculate secondary second-order coefficients
			aNeqSndXY += hermSndXY[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil));
			aNeqSndXZ += hermSndXZ[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil));
			aNeqSndYZ += hermSndYZ[dir] * (getDistribution(alternating, dir) - getEquilibriumDistribution(dir, equil));
		}
	}

	//!set Pi1 tensor elements at node
	Pi1Tensor_[0][0] = aNeqSndXX;
	Pi1Tensor_[1][1] = aNeqSndYY;
	Pi1Tensor_[2][2] = aNeqSndZZ;

	Pi1Tensor_[0][1] = aNeqSndXY;
	Pi1Tensor_[0][2] = aNeqSndXZ;
	Pi1Tensor_[1][2] = aNeqSndYZ;
}


///////////////////////////////////////////////////////////////////
//!MRT	COLLISION OPERATORS
///////////////////////////////////////////////////////////////////

void Node::collideMRT_RMbasis(Control& bc, bool& alternating)
{
	double rawMoments[19] = { 0. };
	this->calcAllMoments_RMbasis(bc, rawMoments, distr_[alternating]); //calc raw moments

	//!set density
	double dens = rawMoments[0];

	double relaxationParameters[19] = { 0. };

	//!set omega to current values
	for (unsigned short dir = 0; dir < 19; dir++) relaxationParameters[dir] = 1.;// this->getOmega();

	//!set omega to current values --> relaxationParameters 5-9: standard MRT; relaxationParameters 4-9=omega -->diverges; relaxationParameters 4-9 & 15-18=omega --> stabil
	//relaxationParameters[4] = this->getCollisionFrequency();
	relaxationParameters[5] = this->getCollisionFrequency();
	relaxationParameters[6] = this->getCollisionFrequency();
	relaxationParameters[7] = this->getCollisionFrequency();
	relaxationParameters[8] = this->getCollisionFrequency();
	relaxationParameters[9] = this->getCollisionFrequency();
	//relaxationParameters[15] = this->getCollisionFrequency();
	//relaxationParameters[16] = this->getCollisionFrequency();
	//relaxationParameters[17] = this->getCollisionFrequency();
	//relaxationParameters[18] = this->getCollisionFrequency();

	//!calc equlibrium
	double eq[19] = { 0. };
	Vector3D velo(rawMoments[1] / rawMoments[0], rawMoments[2] / rawMoments[0], rawMoments[3] / rawMoments[0]);
	this->calcEquilibrium4_RMbasis(eq, rawMoments[0], &velo, bc); //use extended eq

	double eqMoments[19] = { 0. };
	this->calcAllMoments_RMbasis(bc, eqMoments, eq);

	//!do relaxation of moments
	//!moments 0 to 3 (density, momentumX, momentumY, momentumZ) are invariant, i.e. do not change by collision
	for (int dir = 0; dir < 19; dir++) rawMoments[dir] += relaxationParameters[dir] * (eqMoments[dir] - rawMoments[dir]);

	//!do retransformation
	double dummy[19];
	this->transformMomentsToDistributions_RMbasis(dummy, rawMoments, bc); //!calc distributions by re transformation
}

void Node::calcEquilibrium4_RMbasis(double* equil, double density, double* velo, Control& bc)		//!method to calc eq distribution from raw moments
{
	double ux = velo[0];
	double uy = velo[1];
	double uz = velo[2];

	double uxSq = pow(ux, 2);
	double uySq = pow(uy, 2);
	double uzSq = pow(uz, 2);

	double xsi = 1. / bc.getMolecularVelocity();
	double xsiSq = pow(xsi, 2.);
	double xsiPow3 = pow(xsi, 3.);
	double xsiPow4 = pow(xsi, 4.);

	double* weights = bc.getWeightFactors();

	//!COREIXAS' EXTENDED EQUILIBRIUM BASED ON RAW MOMENT FORMALISM
	//!calc eq distributions for all directions, consider higer order extension of equilibrium,
	//!see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. H15 - H21
	equil[0] = density * weights[0] * (1 + 3 * ux * xsi + 3 * (uxSq - uySq - uzSq) * xsiSq - 9 * (ux * uySq + ux * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uxSq * uzSq) * xsiPow4);
	equil[1] = density * weights[1] * (1 + 3 * uy * xsi + 3 * (-uxSq + uySq - uzSq) * xsiSq - 9 * (uy * uxSq + uy * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uySq * uzSq) * xsiPow4);
	equil[2] = density * weights[2] * (1 - 3 * ux * xsi + 3 * (uxSq - uySq - uzSq) * xsiSq + 9 * (ux * uySq + ux * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uxSq * uzSq) * xsiPow4);
	equil[3] = density * weights[3] * (1 - 3 * uy * xsi + 3 * (-uxSq + uySq - uzSq) * xsiSq + 9 * (uy * uxSq + uy * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uySq * uzSq) * xsiPow4);
	equil[4] = density * weights[4] * (1 + 3 * uz * xsi + 3 * (-uxSq - uySq + uzSq) * xsiSq - 9 * (uz * uxSq + uz * uySq) * xsiPow3 - 9 * (uxSq * uzSq + uySq * uzSq) * xsiPow4);
	equil[5] = density * weights[5] * (1 - 3 * uz * xsi + 3 * (-uxSq - uySq + uzSq) * xsiSq + 9 * (uz * uxSq + uz * uySq) * xsiPow3 - 9 * (uxSq * uzSq + uySq * uzSq) * xsiPow4);
	equil[6] = density * weights[6] * (1 + 3 * (ux + uz) * xsi + 3 * (uxSq + uzSq) * xsiSq + 9 * ux * uz * xsiSq + 9 * (uxSq * uz + ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[7] = density * weights[7] * (1 + 3 * (ux - uz) * xsi + 3 * (uxSq + uzSq) * xsiSq - 9 * ux * uz * xsiSq + 9 * (-uxSq * uz + ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[8] = density * weights[8] * (1 + 3 * (uy + uz) * xsi + 3 * (uySq + uzSq) * xsiSq + 9 * uy * uz * xsiSq + 9 * (uySq * uz + uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[9] = density * weights[9] * (1 + 3 * (uy - uz) * xsi + 3 * (uySq + uzSq) * xsiSq - 9 * uy * uz * xsiSq + 9 * (-uySq * uz + uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[10] = density * weights[10] * (1 + 3 * (-ux + uz) * xsi + 3 * (uxSq + uzSq) * xsiSq - 9 * ux * uz * xsiSq + 9 * (uxSq * uz - ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[11] = density * weights[11] * (1 + 3 * (-ux - uz) * xsi + 3 * (uxSq + uzSq) * xsiSq + 9 * ux * uz * xsiSq + 9 * (-uxSq * uz - ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[12] = density * weights[12] * (1 + 3 * (-uy + uz) * xsi + 3 * (uySq + uzSq) * xsiSq - 9 * uy * uz * xsiSq + 9 * (uySq * uz - uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[13] = density * weights[13] * (1 + 3 * (-uy - uz) * xsi + 3 * (uySq + uzSq) * xsiSq + 9 * uy * uz * xsiSq + 9 * (-uySq * uz - uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[14] = density * weights[14] * (1 + 3 * (-ux + uy) * xsi + 3 * (uxSq + uySq) * xsiSq - 9 * ux * uy * xsiSq + 9 * (uxSq * uy - ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[15] = density * weights[15] * (1 + 3 * (ux - uy) * xsi + 3 * (uxSq + uySq) * xsiSq - 9 * ux * uy * xsiSq + 9 * (-uxSq * uy + ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[16] = density * weights[16] * (1 + 3 * (ux + uy) * xsi + 3 * (uxSq + uySq) * xsiSq + 9 * ux * uy * xsiSq + 9 * (uxSq * uy + ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[17] = density * weights[17] * (1 + 3 * (-ux - uy) * xsi + 3 * (uxSq + uySq) * xsiSq + 9 * ux * uy * xsiSq + 9 * (-uxSq * uy - ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[18] = density * weights[18] * (1 - (uxSq + uySq + uzSq) * xsiSq + 3 * (uxSq * uySq + uxSq * uzSq + uySq * uzSq) * xsiPow4);

}//!end calcEquilibrium4_RMbasis


void Node::calcEquilibrium4_RMbasis(double* equil, double density, Vector3D* velo, Control& bc)		//!method to calc eq distribution from raw moments
{
	double ux = velo->getX();
	double uy = velo->getY();
	double uz = velo->getZ();

	double uxSq = pow(ux, 2);
	double uySq = pow(uy, 2);
	double uzSq = pow(uz, 2);

	double xsi = 1. / bc.getMolecularVelocity();
	double xsiSq = pow(xsi, 2.);
	double xsiPow3 = pow(xsi, 3.);
	double xsiPow4 = pow(xsi, 4.);

	double* weights = bc.getWeightFactors();

	//!COREIXAS' EXTENDED EQUILIBRIUM BASED ON RAW MOMENT FORMALISM
	//!calc eq distributions for all directions, consider higer order extension of equilibrium,
	//!see Coreixas et al.: Comprehensive comparison of collision models in the lattice Boltzmann framework: Theoretical investigations. PHYSICAL REVIEW E 100, 033305 (2019), eq. H15 - H21
	equil[0] = density * weights[0] * (1 + 3 * ux * xsi + 3 * (uxSq - uySq - uzSq) * xsiSq - 9 * (ux * uySq + ux * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uxSq * uzSq) * xsiPow4);
	equil[1] = density * weights[1] * (1 + 3 * uy * xsi + 3 * (-uxSq + uySq - uzSq) * xsiSq - 9 * (uy * uxSq + uy * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uySq * uzSq) * xsiPow4);
	equil[2] = density * weights[2] * (1 - 3 * ux * xsi + 3 * (uxSq - uySq - uzSq) * xsiSq + 9 * (ux * uySq + ux * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uxSq * uzSq) * xsiPow4);
	equil[3] = density * weights[3] * (1 - 3 * uy * xsi + 3 * (-uxSq + uySq - uzSq) * xsiSq + 9 * (uy * uxSq + uy * uzSq) * xsiPow3 - 9 * (uxSq * uySq + uySq * uzSq) * xsiPow4);
	equil[4] = density * weights[4] * (1 + 3 * uz * xsi + 3 * (-uxSq - uySq + uzSq) * xsiSq - 9 * (uz * uxSq + uz * uySq) * xsiPow3 - 9 * (uxSq * uzSq + uySq * uzSq) * xsiPow4);
	equil[5] = density * weights[5] * (1 - 3 * uz * xsi + 3 * (-uxSq - uySq + uzSq) * xsiSq + 9 * (uz * uxSq + uz * uySq) * xsiPow3 - 9 * (uxSq * uzSq + uySq * uzSq) * xsiPow4);
	equil[6] = density * weights[6] * (1 + 3 * (ux + uz) * xsi + 3 * (uxSq + uzSq) * xsiSq + 9 * ux * uz * xsiSq + 9 * (uxSq * uz + ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[7] = density * weights[7] * (1 + 3 * (ux - uz) * xsi + 3 * (uxSq + uzSq) * xsiSq - 9 * ux * uz * xsiSq + 9 * (-uxSq * uz + ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[8] = density * weights[8] * (1 + 3 * (uy + uz) * xsi + 3 * (uySq + uzSq) * xsiSq + 9 * uy * uz * xsiSq + 9 * (uySq * uz + uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[9] = density * weights[9] * (1 + 3 * (uy - uz) * xsi + 3 * (uySq + uzSq) * xsiSq - 9 * uy * uz * xsiSq + 9 * (-uySq * uz + uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[10] = density * weights[10] * (1 + 3 * (-ux + uz) * xsi + 3 * (uxSq + uzSq) * xsiSq - 9 * ux * uz * xsiSq + 9 * (uxSq * uz - ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[11] = density * weights[11] * (1 + 3 * (-ux - uz) * xsi + 3 * (uxSq + uzSq) * xsiSq + 9 * ux * uz * xsiSq + 9 * (-uxSq * uz - ux * uzSq) * xsiPow3 + 9 * uxSq * uzSq * xsiPow4);
	equil[12] = density * weights[12] * (1 + 3 * (-uy + uz) * xsi + 3 * (uySq + uzSq) * xsiSq - 9 * uy * uz * xsiSq + 9 * (uySq * uz - uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[13] = density * weights[13] * (1 + 3 * (-uy - uz) * xsi + 3 * (uySq + uzSq) * xsiSq + 9 * uy * uz * xsiSq + 9 * (-uySq * uz - uy * uzSq) * xsiPow3 + 9 * uySq * uzSq * xsiPow4);
	equil[14] = density * weights[14] * (1 + 3 * (-ux + uy) * xsi + 3 * (uxSq + uySq) * xsiSq - 9 * ux * uy * xsiSq + 9 * (uxSq * uy - ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[15] = density * weights[15] * (1 + 3 * (ux - uy) * xsi + 3 * (uxSq + uySq) * xsiSq - 9 * ux * uy * xsiSq + 9 * (-uxSq * uy + ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[16] = density * weights[16] * (1 + 3 * (ux + uy) * xsi + 3 * (uxSq + uySq) * xsiSq + 9 * ux * uy * xsiSq + 9 * (uxSq * uy + ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[17] = density * weights[17] * (1 + 3 * (-ux - uy) * xsi + 3 * (uxSq + uySq) * xsiSq + 9 * ux * uy * xsiSq + 9 * (-uxSq * uy - ux * uySq) * xsiPow3 + 9 * uxSq * uySq * xsiPow4);
	equil[18] = density * weights[18] * (1 - (uxSq + uySq + uzSq) * xsiSq + 3 * (uxSq * uySq + uxSq * uzSq + uySq * uzSq) * xsiPow4);

}//!end calcEquilibrium4_RMbasis


void Node::calcAllMoments_RMbasis(Control& bc, double* moments, double* distribution) const   //!choose distribution by alternating (double (*functionParameters)[19],double* moments,bool alternating) const
{
	double molecularVelocity = bc.getMolecularVelocity();
	double molecularVelocitySquared = pow(molecularVelocity, 2.0);//bc.getStatics()->getMolecularVelocitySquared();
	double molecularVelocityPower3 = pow(molecularVelocity, 3.0);//bc.getStatics()->getMolecularVelocityPower3();
	double molecularVelocityPower4 = pow(molecularVelocity, 4.0);//bc.getStatics()->getMolecularVelocityPower4();

    //!calc moments --> save operations by hard coding
	moments[0] = 0.;
	for (unsigned short direction = 0; direction < 19; direction++) moments[0] += distribution[direction];	//!density=M_000

	moments[1] = distribution[0] - distribution[2] + distribution[6] + distribution[7] - distribution[10] - distribution[11] - distribution[14] + distribution[15] + distribution[16] - distribution[17];
	moments[1] *= molecularVelocity;	//!momentum in x-direction=M_100

	moments[2] = distribution[1] - distribution[3] + distribution[8] + distribution[9] - distribution[12] - distribution[13] + distribution[14] - distribution[15] + distribution[16] - distribution[17];
	moments[2] *= molecularVelocity;	//!momentum in y-direction=M_010

	moments[3] = distribution[4] - distribution[5] + distribution[6] - distribution[7] + distribution[8] - distribution[9] + distribution[10] - distribution[11] + distribution[12] - distribution[13];
	moments[3] *= molecularVelocity;	//!momentum in z-direction=M_001

	moments[4] = 0;
	for (unsigned short direction = 0; direction < 6; direction++) moments[4] += distribution[direction];
	double temp = 0.;
	for (unsigned short direction = 6; direction < 18; direction++) temp += distribution[direction];
	moments[4] += 2. * temp;
	moments[4] *= molecularVelocitySquared;	//!trace T=M_200+M_020+M_002 of stress tensor, i.e. Pi_xx+Pi_yy+Pi_zz (in lattice units)

	moments[5] = distribution[0] + distribution[2] - distribution[4] - distribution[5] - distribution[8] - distribution[9] - distribution[12] - distribution[13] + distribution[14] + distribution[15] + distribution[16] + distribution[17];
	moments[5] *= molecularVelocitySquared;	//!N_xz=M_200-M_002, i.e. Pi_xx- Pi_zz (in lattice units)

	moments[6] = distribution[1] + distribution[3] - distribution[4] - distribution[5] - distribution[6] - distribution[7] - distribution[10] - distribution[11] + distribution[14] + distribution[15] + distribution[16] + distribution[17];
	moments[6] *= molecularVelocitySquared;	//!N_yz=M_020-M_002, i.e. Pi_yy- Pi_zz (in lattice units)

	moments[7] = -distribution[14] - distribution[15] + distribution[16] + distribution[17];
	moments[7] *= molecularVelocitySquared;	//!shear Pi_xy=M_110 (in lattice units)

	moments[8] = distribution[6] - distribution[7] - distribution[10] + distribution[11];
	moments[8] *= molecularVelocitySquared;	//!shear Pi_xz=M_101 (in lattice units)

	moments[9] = distribution[8] - distribution[9] - distribution[12] + distribution[13];
	moments[9] *= molecularVelocitySquared;	//!shear Pi_yz=M_011 (in lattice units)

	moments[10] = distribution[6] - distribution[7] + distribution[10] - distribution[11];
	moments[10] *= molecularVelocityPower3;	//!heat flux Q_xxz=M_201 (in lattice units)

	moments[11] = distribution[6] + distribution[7] - distribution[10] - distribution[11];
	moments[11] *= molecularVelocityPower3;	//!heat flux Q_xzz=M_102 (in lattice units)

	moments[12] = distribution[8] - distribution[9] + distribution[12] - distribution[13];
	moments[12] *= molecularVelocityPower3;	//!heat flux Q_yyz=M_021 (in lattice units)

	moments[13] = distribution[8] + distribution[9] - distribution[12] - distribution[13];
	moments[13] *= molecularVelocityPower3;	//!heat flux Q_yzz=M_012 (in lattice units)

	moments[14] = distribution[14] - distribution[15] + distribution[16] - distribution[17];
	moments[14] *= molecularVelocityPower3;	//!heat flux Q_xxy=M_210 (in lattice units)

	moments[15] = -distribution[14] + distribution[15] + distribution[16] - distribution[17];
	moments[15] *= molecularVelocityPower3;	//!heat flux Q_xyy=M_120 (in lattice units)

	moments[16] = distribution[6] + distribution[7] + distribution[10] + distribution[11];
	moments[16] *= molecularVelocityPower4;	//!M_202 (in lattice units)

	moments[17] = distribution[14] + distribution[15] + distribution[16] + distribution[17];
	moments[17] *= molecularVelocityPower4;	//!M_220 (in lattice units)

	moments[18] = distribution[8] + distribution[9] + distribution[12] + distribution[13];
	moments[18] *= molecularVelocityPower4;	//!M_022 (in lattice units)

}//!end calcAllMoments_RMbasis


void Node::transformMomentsToDistributions_RMbasis(double* distribution, double* moments, Control& bc)		//!method to transform raw moments to distributions
{
	//!inverse molecular velocities
	static double xsi = 1. / bc.getMolecularVelocity();
	static double xsiPower2 = pow(xsi, 2);
	static double xsiPower3 = pow(xsi, 3);
	static double xsiPower4 = pow(xsi, 4);

	static double oneDivSix = 1. / 6.;

	//!transform moments to distributions for all directions (do only necessary calculations, i.e. matrix entry !=0)
	distribution[0] = 0.5 * moments[1] * xsi + oneDivSix * (moments[4] + 2. * moments[5] - moments[6]) * xsiPower2 - 0.5 * (moments[11] + moments[15]) * xsiPower3 - 0.5 * (moments[16] + moments[17]) * xsiPower4;
	distribution[1] = 0.5 * moments[2] * xsi + oneDivSix * (moments[4] - moments[5] + 2. * moments[6]) * xsiPower2 - 0.5 * (moments[13] + moments[14]) * xsiPower3 - 0.5 * (moments[17] + moments[18]) * xsiPower4;
	distribution[2] = -0.5 * moments[1] * xsi + oneDivSix * (moments[4] + 2. * moments[5] - moments[6]) * xsiPower2 + 0.5 * (moments[11] + moments[15]) * xsiPower3 - 0.5 * (moments[16] + moments[17]) * xsiPower4;
	distribution[3] = -0.5 * moments[2] * xsi + oneDivSix * (moments[4] - moments[5] + 2. * moments[6]) * xsiPower2 + 0.5 * (moments[13] + moments[14]) * xsiPower3 - 0.5 * (moments[17] + moments[18]) * xsiPower4;
	distribution[4] = 0.5 * moments[3] * xsi + oneDivSix * (moments[4] - moments[5] - moments[6]) * xsiPower2 - 0.5 * (moments[10] + moments[12]) * xsiPower3 - 0.5 * (moments[16] + moments[18]) * xsiPower4;
	distribution[5] = -0.5 * moments[3] * xsi + oneDivSix * (moments[4] - moments[5] - moments[6]) * xsiPower2 + 0.5 * (moments[10] + moments[12]) * xsiPower3 - 0.5 * (moments[16] + moments[18]) * xsiPower4;
	distribution[6] = 0.25 * moments[8] * xsiPower2 + 0.25 * (moments[10] + moments[11]) * xsiPower3 + 0.25 * moments[16] * xsiPower4;
	distribution[7] = -0.25 * moments[8] * xsiPower2 - 0.25 * (moments[10] - moments[11]) * xsiPower3 + 0.25 * moments[16] * xsiPower4;
	distribution[8] = 0.25 * moments[9] * xsiPower2 + 0.25 * (moments[12] + moments[13]) * xsiPower3 + 0.25 * moments[18] * xsiPower4;
	distribution[9] = -0.25 * moments[9] * xsiPower2 - 0.25 * (moments[12] - moments[13]) * xsiPower3 + 0.25 * moments[18] * xsiPower4;
	distribution[10] = -0.25 * moments[8] * xsiPower2 + 0.25 * (moments[10] - moments[11]) * xsiPower3 + 0.25 * moments[16] * xsiPower4;
	distribution[11] = 0.25 * moments[8] * xsiPower2 - 0.25 * (moments[10] + moments[11]) * xsiPower3 + 0.25 * moments[16] * xsiPower4;
	distribution[12] = -0.25 * moments[9] * xsiPower2 + 0.25 * (moments[12] - moments[13]) * xsiPower3 + 0.25 * moments[18] * xsiPower4;
	distribution[13] = 0.25 * moments[9] * xsiPower2 - 0.25 * (moments[12] + moments[13]) * xsiPower3 + 0.25 * moments[18] * xsiPower4;
	distribution[14] = -0.25 * moments[7] * xsiPower2 + 0.25 * (moments[14] - moments[15]) * xsiPower3 + 0.25 * moments[17] * xsiPower4;
	distribution[15] = -0.25 * moments[7] * xsiPower2 - 0.25 * (moments[14] - moments[15]) * xsiPower3 + 0.25 * moments[17] * xsiPower4;
	distribution[16] = 0.25 * moments[7] * xsiPower2 + 0.25 * (moments[14] + moments[15]) * xsiPower3 + 0.25 * moments[17] * xsiPower4;
	distribution[17] = 0.25 * moments[7] * xsiPower2 - 0.25 * (moments[14] + moments[15]) * xsiPower3 + 0.25 * moments[17] * xsiPower4;
	distribution[18] = moments[0] - moments[4] * xsiPower2 + (moments[16] + moments[17] + moments[18]) * xsiPower4;
}//!end transformRawMomentsToDistributions


void Node::collideMRT_GSbasis(Control& bc, bool& alternating)
{
	int momentsToRelax[15] = { 1,2,4,6,8,9,10,11,12,13,14,15,16,17,18 };
	int momentsNotToRelax[4] = { 0,3,5,7 };
	//!enter values for relaxation parameters
	double relaxationParameters[19] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };

	//!calc moments at this time step
	double moments[19];	//!arrangement of moment is the same as basis vectors (i.e (density,energy,energy squared,momentumX,heatFluxX,momentumY,heatFluxY,momentumZ,heatFluxZ,pressureXX,shearStressXX,pressureZZ,shearStressZZ,pressureXY,pressureYZ,pressureXZ,mx,my,mz)
	this->calcAllMoments_GSbasis(bc, moments, alternating);

	//!calc equilibrium moments
	double equilibriumMoments[19] = { 0. };
	this->calcEquilibriumMoments2_GSbasis(equilibriumMoments, moments[0], moments[3], moments[5], moments[7]);

	//!collision invariants are not affected
	relaxationParameters[0] = 0.;
	relaxationParameters[3] = 0.;
	relaxationParameters[5] = 0.;
	relaxationParameters[7] = 0.;

	//!set omega to current values (necessary for multi level & turbulent flows)
	relaxationParameters[1] = this->getCollisionFrequency(); //!set bulk viscosity equal to shear viscosity for aeracoustic simulations
	relaxationParameters[9] = relaxationParameters[1];
	relaxationParameters[11] = relaxationParameters[1];
	relaxationParameters[13] = relaxationParameters[1];
	relaxationParameters[14] = relaxationParameters[1];
	relaxationParameters[15] = relaxationParameters[1];

	//!do relaxation for necessary moments
	double changeOfMoments[19];
	for (int dir = 0; dir < 15; dir++) changeOfMoments[momentsToRelax[dir]] = relaxationParameters[momentsToRelax[dir]] * (equilibriumMoments[momentsToRelax[dir]] - moments[momentsToRelax[dir]]);
	for (int dir = 0; dir < 4; dir++)  changeOfMoments[momentsNotToRelax[dir]] = 0.;

	//!do retransformation
	double dummy[19];
	this->transformMomentsToDistributions_GSbasis(dummy, changeOfMoments, bc);
}


void Node::calcEquilibriumMoments2_GSbasis(double* equilibriumMoments, double& density, double& momentumX, double& momentumY, double& momentumZ)
{
	//!calc equilibrium moments
	//!set conserved moments
	equilibriumMoments[0] = density;		//!density
	equilibriumMoments[3] = momentumX;		//!momentum x
	equilibriumMoments[5] = momentumY;		//!momentum y
	equilibriumMoments[7] = momentumZ;		//!momentum z

	//!set conserved moments --> all other equilibrium moments are 0
	double reverseDensity = 1. / equilibriumMoments[0];
	double momentumXsquared = equilibriumMoments[3] * equilibriumMoments[3];
	double momentumYsquared = equilibriumMoments[5] * equilibriumMoments[5];
	double momentumZsquared = equilibriumMoments[7] * equilibriumMoments[7];
	equilibriumMoments[1] = (momentumXsquared + momentumYsquared + momentumZsquared) * reverseDensity;
	equilibriumMoments[2] = 0.;
	equilibriumMoments[4] = 0.;
	equilibriumMoments[6] = 0.;
	equilibriumMoments[8] = 0.;
	equilibriumMoments[9] = (2. * momentumXsquared - momentumYsquared - momentumZsquared) * reverseDensity;
	equilibriumMoments[10] = 0.;
	equilibriumMoments[11] = (momentumYsquared - momentumZsquared) * reverseDensity;
	equilibriumMoments[12] = 0.;
	equilibriumMoments[13] = equilibriumMoments[3] * equilibriumMoments[5] * reverseDensity;
	equilibriumMoments[14] = equilibriumMoments[5] * equilibriumMoments[7] * reverseDensity;
	equilibriumMoments[15] = equilibriumMoments[3] * equilibriumMoments[7] * reverseDensity;
	equilibriumMoments[16] = 0.;
	equilibriumMoments[17] = 0.;
	equilibriumMoments[18] = 0.;
}


void Node::transformMomentsToDistributions_GSbasis(double* distributions, double* moments, Control& bc)
{
	double xsi = 1. / bc.getMolecularVelocity();
	double xsiPower2 = pow(xsi, 2.);
	double xsiPower3 = pow(xsi, 3.);
	double xsiPower4 = pow(xsi, 4.);

	//!multiplicators for mrt back transformation
	double prefactors[19][19] = { {1. / 18.,1. / 18.,1. / 18.,1. / 18.,1. / 18.,1. / 18.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 36.,1. / 3.},
									  {0.,0.,0.,0.,0.,0.,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,1. / 24. * xsiPower2,-0.5 * xsiPower2},
									  {-1. / 18. * xsiPower4,-1. / 18. * xsiPower4,-1. / 18. * xsiPower4,-1. / 18. * xsiPower4,-1. / 18. * xsiPower4,-1. / 18. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 72. * xsiPower4,1. / 6. * xsiPower4},
									  {1. / 6. * xsi,0.,-1. / 6. * xsi,0.,0.,0.,1. / 12. * xsi,1. / 12. * xsi,0.,0.,-1. / 12. * xsi,-1. / 12. * xsi,0.,0.,-1. / 12. * xsi,1. / 12. * xsi,1. / 12. * xsi,-1. / 12. * xsi,0.},
									  {-1. / 6. * xsiPower3,0.,1. / 6. * xsiPower3,0.,0.,0.,1. / 24. * xsiPower3,1. / 24. * xsiPower3,0.,0.,-1. / 24. * xsiPower3,-1. / 24. * xsiPower3,0.,0.,-1. / 24. * xsiPower3,1. / 24. * xsiPower3,1. / 24. * xsiPower3,-1. / 24. * xsiPower3,0.},
									  {0.,1. / 6. * xsi,0.,-1. / 6. * xsi,0.,0.,0.,0.,1. / 12. * xsi,1. / 12. * xsi,0.,0.,-1. / 12. * xsi,-1. / 12. * xsi,1. / 12. * xsi,-1. / 12. * xsi,1. / 12. * xsi,-1. / 12. * xsi,0.},
									  {0.,-1. / 6. * xsiPower3,0.,1. / 6. * xsiPower3,0.,0.,0.,0.,1. / 24. * xsiPower3,1. / 24. * xsiPower3,0.,0.,-1. / 24. * xsiPower3,-1. / 24. * xsiPower3,1. / 24. * xsiPower3,-1. / 24. * xsiPower3,1. / 24. * xsiPower3,-1. / 24. * xsiPower3,0.},
									  {0.,0.,0.,0.,1. / 6. * xsi,-1. / 6. * xsi,1. / 12. * xsi,-1. / 12. * xsi,1. / 12. * xsi,-1. / 12. * xsi,1. / 12. * xsi,-1. / 12. * xsi,1. / 12. * xsi,-1. / 12. * xsi,0.,0.,0.,0.,0.},
									  {0.,0.,0.,0.,-1. / 6. * xsiPower3,1. / 6. * xsiPower3,1. / 24. * xsiPower3,-1. / 24. * xsiPower3,1. / 24. * xsiPower3,-1. / 24. * xsiPower3,1. / 24. * xsiPower3,-1. / 24. * xsiPower3,1. / 24. * xsiPower3,-1. / 24. * xsiPower3,0.,0.,0.,0.,0.},
									  {1. / 12. * xsiPower2,-1. / 24. * xsiPower2,1. / 12. * xsiPower2,-1. / 24. * xsiPower2,-1. / 24. * xsiPower2,-1. / 24. * xsiPower2,1. / 48. * xsiPower2,1. / 48. * xsiPower2,-1. / 24. * xsiPower2,-1. / 24. * xsiPower2,1. / 48. * xsiPower2,1. / 48. * xsiPower2,-1. / 24. * xsiPower2,-1. / 24. * xsiPower2,1. / 48. * xsiPower2,1. / 48. * xsiPower2,1. / 48. * xsiPower2,1. / 48. * xsiPower2,0.},
									  {-1. / 12. * xsiPower4,1. / 24. * xsiPower4,-1. / 12. * xsiPower4,1. / 24. * xsiPower4,1. / 24. * xsiPower4,1. / 24. * xsiPower4,1. / 48. * xsiPower4,1. / 48. * xsiPower4,-1. / 24. * xsiPower4,-1. / 24. * xsiPower4,1. / 48. * xsiPower4,1. / 48. * xsiPower4,-1. / 24. * xsiPower4,-1. / 24. * xsiPower4,1. / 48. * xsiPower4,1. / 48. * xsiPower4,1. / 48. * xsiPower4,1. / 48. * xsiPower4,0.},
									  {0.,1. / 8. * xsiPower2,0.,1. / 8. * xsiPower2,-1. / 8. * xsiPower2,-1. / 8. * xsiPower2,-1. / 16. * xsiPower2,-1. / 16. * xsiPower2,0.,0.,-1. / 16. * xsiPower2,-1. / 16. * xsiPower2,0.,0.,1. / 16. * xsiPower2,1. / 16. * xsiPower2,1. / 16. * xsiPower2,1. / 16. * xsiPower2,0.},
									  {0.,-1. / 8. * xsiPower4,0.,-1. / 8. * xsiPower4,1. / 8. * xsiPower4,1. / 8. * xsiPower4,-1. / 16. * xsiPower4,-1. / 16. * xsiPower4,0.,0.,-1. / 16. * xsiPower4,-1. / 16. * xsiPower4,0.,0.,1. / 16. * xsiPower4,1. / 16. * xsiPower4,1. / 16. * xsiPower4,1. / 16. * xsiPower4,0.},
									  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1. / 4 * xsiPower2,-1. / 4 * xsiPower2,1. / 4 * xsiPower2,1. / 4 * xsiPower2,0.},
									  {0.,0.,0.,0.,0.,0.,0.,0.,1. / 4 * xsiPower2,-1. / 4 * xsiPower2,0.,0.,-1. / 4 * xsiPower2,1. / 4 * xsiPower2,0.,0.,0.,0.,0.},
									  {0.,0.,0.,0.,0.,0.,1. / 4 * xsiPower2,-1. / 4 * xsiPower2,0.,0.,-1. / 4 * xsiPower2,1. / 4 * xsiPower2,0.,0.,0.,0.,0.,0.,0.},
									  {0.,0.,0.,0.,0.,0.,-1. / 8. * xsiPower3,-1. / 8. * xsiPower3,0.,0.,1. / 8. * xsiPower3,1. / 8. * xsiPower3,0.,0.,-1. / 8. * xsiPower3,1. / 8. * xsiPower3,1. / 8. * xsiPower3,-1. / 8. * xsiPower3,0.},
									  {0.,0.,0.,0.,0.,0.,0.,0.,1. / 8. * xsiPower3,1. / 8. * xsiPower3,0.,0.,-1. / 8. * xsiPower3,-1. / 8. * xsiPower3,-1. / 8. * xsiPower3,1. / 8. * xsiPower3,-1. / 8. * xsiPower3,1. / 8. * xsiPower3,0.},
									  {0.,0.,0.,0.,0.,0.,1. / 8. * xsiPower3,-1. / 8. * xsiPower3,-1. / 8. * xsiPower3,1. / 8. * xsiPower3,1. / 8. * xsiPower3,-1. / 8. * xsiPower3,-1. / 8. * xsiPower3,1. / 8. * xsiPower3,0.,0.,0.,0.,0.}
	};

	//!transform moments to distributions for all directions (do only necessary calculations, i.e. matrix entry !=0)
	distributions[0] = prefactors[0][0] * moments[0] + prefactors[2][0] * moments[2] + prefactors[3][0] * moments[3] + prefactors[4][0] * moments[4] + prefactors[9][0] * moments[9] + prefactors[10][0] * moments[10];
	distributions[1] = prefactors[0][1] * moments[0] + prefactors[2][1] * moments[2] + prefactors[5][1] * moments[5] + prefactors[6][1] * moments[6] + prefactors[9][1] * moments[9] + prefactors[10][1] * moments[10] + prefactors[11][1] * moments[11] + prefactors[12][1] * moments[12];
	distributions[2] = prefactors[0][2] * moments[0] + prefactors[2][2] * moments[2] + prefactors[3][2] * moments[3] + prefactors[4][2] * moments[4] + prefactors[9][2] * moments[9] + prefactors[10][2] * moments[10];
	distributions[3] = prefactors[0][3] * moments[0] + prefactors[2][3] * moments[2] + prefactors[5][3] * moments[5] + prefactors[6][3] * moments[6] + prefactors[9][3] * moments[9] + prefactors[10][3] * moments[10] + prefactors[11][3] * moments[11] + prefactors[12][3] * moments[12];
	distributions[4] = prefactors[0][4] * moments[0] + prefactors[2][4] * moments[2] + prefactors[7][4] * moments[7] + prefactors[8][4] * moments[8] + prefactors[9][4] * moments[9] + prefactors[10][4] * moments[10] + prefactors[11][4] * moments[11] + prefactors[12][4] * moments[12];
	distributions[5] = prefactors[0][5] * moments[0] + prefactors[2][5] * moments[2] + prefactors[7][5] * moments[7] + prefactors[8][5] * moments[8] + prefactors[9][5] * moments[9] + prefactors[10][5] * moments[10] + prefactors[11][5] * moments[11] + prefactors[12][5] * moments[12];
	distributions[6] = prefactors[0][6] * moments[0] + prefactors[1][6] * moments[1] + prefactors[2][6] * moments[2] + prefactors[3][6] * moments[3] + prefactors[4][6] * moments[4] + prefactors[7][6] * moments[7] + prefactors[8][6] * moments[8] + prefactors[9][6] * moments[9] + prefactors[10][6] * moments[10]
		+ prefactors[11][6] * moments[11] + prefactors[12][6] * moments[12] + prefactors[15][6] * moments[15] + prefactors[16][6] * moments[16] + prefactors[18][6] * moments[18];
	distributions[7] = prefactors[0][7] * moments[0] + prefactors[1][7] * moments[1] + prefactors[2][7] * moments[2] + prefactors[3][7] * moments[3] + prefactors[4][7] * moments[4] + prefactors[7][7] * moments[7] + prefactors[8][7] * moments[8] + prefactors[9][7] * moments[9] + prefactors[10][7] * moments[10]
		+ prefactors[11][7] * moments[11] + prefactors[12][7] * moments[12] + prefactors[15][7] * moments[15] + prefactors[16][7] * moments[16] + prefactors[18][7] * moments[18];
	distributions[8] = prefactors[0][8] * moments[0] + prefactors[1][8] * moments[1] + prefactors[2][8] * moments[2] + prefactors[5][8] * moments[5] + prefactors[6][8] * moments[6] + prefactors[7][8] * moments[7] + prefactors[8][8] * moments[8] + prefactors[9][8] * moments[9] + prefactors[10][8] * moments[10]
		+ prefactors[14][8] * moments[14] + prefactors[17][8] * moments[17] + prefactors[18][8] * moments[18];
	distributions[9] = prefactors[0][9] * moments[0] + prefactors[1][9] * moments[1] + prefactors[2][9] * moments[2] + prefactors[5][9] * moments[5] + prefactors[6][9] * moments[6] + prefactors[7][9] * moments[7] + prefactors[8][9] * moments[8] + prefactors[9][9] * moments[9] + prefactors[10][9] * moments[10]
		+ prefactors[14][9] * moments[14] + prefactors[17][9] * moments[17] + prefactors[18][9] * moments[18];
	distributions[10] = prefactors[0][10] * moments[0] + prefactors[1][10] * moments[1] + prefactors[2][10] * moments[2] + prefactors[3][10] * moments[3] + prefactors[4][10] * moments[4] + prefactors[7][10] * moments[7] + prefactors[8][10] * moments[8] + prefactors[9][10] * moments[9] + prefactors[10][10] * moments[10]
		+ prefactors[11][10] * moments[11] + prefactors[12][10] * moments[12] + prefactors[15][10] * moments[15] + prefactors[16][10] * moments[16] + prefactors[18][10] * moments[18];
	distributions[11] = prefactors[0][11] * moments[0] + prefactors[1][11] * moments[1] + prefactors[2][11] * moments[2] + prefactors[3][11] * moments[3] + prefactors[4][11] * moments[4] + prefactors[7][11] * moments[7] + prefactors[8][11] * moments[8] + prefactors[9][11] * moments[9] + prefactors[10][11] * moments[10]
		+ prefactors[11][11] * moments[11] + prefactors[12][11] * moments[12] + prefactors[15][11] * moments[15] + prefactors[16][11] * moments[16] + prefactors[18][11] * moments[18];
	distributions[12] = prefactors[0][12] * moments[0] + prefactors[1][12] * moments[1] + prefactors[2][12] * moments[2] + prefactors[5][12] * moments[5] + prefactors[6][12] * moments[6] + prefactors[7][12] * moments[7] + prefactors[8][12] * moments[8] + prefactors[9][12] * moments[9] + prefactors[10][12] * moments[10]
		+ prefactors[14][12] * moments[14] + prefactors[17][12] * moments[17] + prefactors[18][12] * moments[18];
	distributions[13] = prefactors[0][13] * moments[0] + prefactors[1][13] * moments[1] + prefactors[2][13] * moments[2] + prefactors[5][13] * moments[5] + prefactors[6][13] * moments[6] + prefactors[7][13] * moments[7] + prefactors[8][13] * moments[8] + prefactors[9][13] * moments[9] + prefactors[10][13] * moments[10]
		+ prefactors[14][13] * moments[14] + prefactors[17][13] * moments[17] + prefactors[18][13] * moments[18];
	distributions[14] = prefactors[0][14] * moments[0] + prefactors[1][14] * moments[1] + prefactors[2][14] * moments[2] + prefactors[3][14] * moments[3] + prefactors[4][14] * moments[4] + prefactors[5][14] * moments[5] + prefactors[6][14] * moments[6] + prefactors[9][14] * moments[9] + prefactors[10][14] * moments[10] + prefactors[11][14] * moments[11]
		+ prefactors[12][14] * moments[12] + prefactors[13][14] * moments[13] + prefactors[16][14] * moments[16] + prefactors[17][14] * moments[17];
	distributions[15] = prefactors[0][15] * moments[0] + prefactors[1][15] * moments[1] + prefactors[2][15] * moments[2] + prefactors[3][15] * moments[3] + prefactors[4][15] * moments[4] + prefactors[5][15] * moments[5] + prefactors[6][15] * moments[6] + prefactors[9][15] * moments[9] + prefactors[10][15] * moments[10] + prefactors[11][15] * moments[11]
		+ prefactors[12][15] * moments[12] + prefactors[13][15] * moments[13] + prefactors[16][15] * moments[16] + prefactors[17][15] * moments[17];
	distributions[16] = prefactors[0][16] * moments[0] + prefactors[1][16] * moments[1] + prefactors[2][16] * moments[2] + prefactors[3][16] * moments[3] + prefactors[4][16] * moments[4] + prefactors[5][16] * moments[5] + prefactors[6][16] * moments[6] + prefactors[9][16] * moments[9] + prefactors[10][16] * moments[10] + prefactors[11][16] * moments[11]
		+ prefactors[12][16] * moments[12] + prefactors[13][16] * moments[13] + prefactors[16][16] * moments[16] + prefactors[17][16] * moments[17];
	distributions[17] = prefactors[0][17] * moments[0] + prefactors[1][17] * moments[1] + prefactors[2][17] * moments[2] + prefactors[3][17] * moments[3] + prefactors[4][17] * moments[4] + prefactors[5][17] * moments[5] + prefactors[6][17] * moments[6] + prefactors[9][17] * moments[9] + prefactors[10][17] * moments[10] + prefactors[11][17] * moments[11]
		+ prefactors[12][17] * moments[12] + prefactors[13][17] * moments[13] + prefactors[16][17] * moments[16] + prefactors[17][17] * moments[17];
	distributions[18] = prefactors[0][18] * moments[0] + prefactors[1][18] * moments[1] + prefactors[2][18] * moments[2];
}


void Node::calcAllMoments_GSbasis(Control& bc, double* moments, bool& alternating)
{
	double molecularVelocity = bc.getMolecularVelocity();
	double molecularVelocitySquared = pow(molecularVelocity, 2.);
	double molecularVelocityCubic = pow(molecularVelocity, 3.);
	double molecularVelocityQuartic = pow(molecularVelocity, 4.);

	//!calc moments --> save operations by hard coding
	moments[0] = 0.;
	for (int dir = 0; dir < 19; dir++) moments[0] += distr_[alternating][dir];	//!density

	moments[1] = 0.;
	for (int dir = 6; dir < 18; dir++) moments[1] += distr_[alternating][dir];
	moments[1] -= distr_[alternating][18];
	moments[1] *= molecularVelocitySquared;	//!function for developing energy

	moments[2] = 0.;
	for (int dir = 0; dir < 6; dir++) moments[2] -= 2. * distr_[alternating][dir];
	for (int dir = 6; dir < 19; dir++) moments[2] += distr_[alternating][dir];
	moments[2] *= molecularVelocityQuartic;	//function for developing energy power2

	moments[3] = distr_[alternating][0] - distr_[alternating][2] + distr_[alternating][6] + distr_[alternating][7] - distr_[alternating][10] - distr_[alternating][11] - distr_[alternating][14] + distr_[alternating][15] + distr_[alternating][16] - distr_[alternating][17];
	moments[3] *= molecularVelocity;	//!function for developing momentum in x-direction

	moments[4] = -2. * distr_[alternating][0] + 2 * distr_[alternating][2] + distr_[alternating][6] + distr_[alternating][7] - distr_[alternating][10] - distr_[alternating][11] - distr_[alternating][14] + distr_[alternating][15] + distr_[alternating][16] - distr_[alternating][17];
	moments[4] *= molecularVelocityCubic;	//!function for developing heat flux in x-direction

	moments[5] = distr_[alternating][1] - distr_[alternating][3] + distr_[alternating][8] + distr_[alternating][9] - distr_[alternating][12] - distr_[alternating][13] + distr_[alternating][14] - distr_[alternating][15] + distr_[alternating][16] - distr_[alternating][17];
	moments[5] *= molecularVelocity;	//!function for developing momentum in y-direction

	moments[6] = -2. * distr_[alternating][1] + 2 * distr_[alternating][3] + distr_[alternating][8] + distr_[alternating][9] - distr_[alternating][12] - distr_[alternating][13] + distr_[alternating][14] - distr_[alternating][15] + distr_[alternating][16] - distr_[alternating][17];
	moments[6] *= molecularVelocityCubic;	//!function for developing heat flux in y-direction

	moments[7] = distr_[alternating][4] - distr_[alternating][5] + distr_[alternating][6] - distr_[alternating][7] + distr_[alternating][8] - distr_[alternating][9] + distr_[alternating][10] - distr_[alternating][11] + distr_[alternating][12] - distr_[alternating][13];
	moments[7] *= molecularVelocity;	//!function for developing momentum in z-direction

	moments[8] = -2. * distr_[alternating][4] + 2 * distr_[alternating][5] + distr_[alternating][6] - distr_[alternating][7] + distr_[alternating][8] - distr_[alternating][9] + distr_[alternating][10] - distr_[alternating][11] + distr_[alternating][12] - distr_[alternating][13];
	moments[8] *= molecularVelocityCubic;	//!function for developing heat flux in z-direction

	moments[9] = 2. * distr_[alternating][0] - distr_[alternating][1] + 2. * distr_[alternating][2] - distr_[alternating][3] - distr_[alternating][4] - distr_[alternating][5] + distr_[alternating][6] + distr_[alternating][7] - 2. * distr_[alternating][8] - 2. * distr_[alternating][9] + distr_[alternating][10] + distr_[alternating][11]
		- 2. * distr_[alternating][12] - 2. * distr_[alternating][13] + distr_[alternating][14] + distr_[alternating][15] + distr_[alternating][16] + distr_[alternating][17];
	moments[9] *= molecularVelocitySquared;	//!function for developing pressure xx

	moments[10] = -2. * distr_[alternating][0] + distr_[alternating][1] - 2. * distr_[alternating][2] + distr_[alternating][3] + distr_[alternating][4] + distr_[alternating][5] + distr_[alternating][6] + distr_[alternating][7] - 2. * distr_[alternating][8] - 2. * distr_[alternating][9] + distr_[alternating][10] + distr_[alternating][11]
		- 2. * distr_[alternating][12] - 2. * distr_[alternating][13] + distr_[alternating][14] + distr_[alternating][15] + distr_[alternating][16] + distr_[alternating][17];
	moments[10] *= molecularVelocityQuartic;	//function for developing shear stress xx

	moments[11] = distr_[alternating][1] + distr_[alternating][3] - distr_[alternating][4] - distr_[alternating][5] - distr_[alternating][6] - distr_[alternating][7] - distr_[alternating][10] - distr_[alternating][11] + distr_[alternating][14] + distr_[alternating][15] + distr_[alternating][16] + distr_[alternating][17];
	moments[11] *= molecularVelocitySquared;	//!function for developing pressure zz

	moments[12] = -distr_[alternating][1] - distr_[alternating][3] + distr_[alternating][4] + distr_[alternating][5] - distr_[alternating][6] - distr_[alternating][7] - distr_[alternating][10] - distr_[alternating][11] + distr_[alternating][14] + distr_[alternating][15] + distr_[alternating][16] + distr_[alternating][17];
	moments[12] *= molecularVelocityQuartic;	//!function for developing shear stress zz

	moments[13] = -distr_[alternating][14] - distr_[alternating][15] + distr_[alternating][16] + distr_[alternating][17];
	moments[13] *= molecularVelocitySquared; //!function for developing pressure xy

	moments[14] = distr_[alternating][8] - distr_[alternating][9] - distr_[alternating][12] + distr_[alternating][13];
	moments[14] *= molecularVelocitySquared; //!function for developing pressure yz

	moments[15] = distr_[alternating][6] - distr_[alternating][7] - distr_[alternating][10] + distr_[alternating][11];
	moments[15] *= molecularVelocitySquared; //!function for developing pressure xz

	moments[16] = -distr_[alternating][6] - distr_[alternating][7] + distr_[alternating][10] + distr_[alternating][11] - distr_[alternating][14] + distr_[alternating][15] + distr_[alternating][16] - distr_[alternating][17];
	moments[16] *= molecularVelocityCubic;	//!function for developing moment mx

	moments[17] = distr_[alternating][8] + distr_[alternating][9] - distr_[alternating][12] - distr_[alternating][13] - distr_[alternating][14] + distr_[alternating][15] - distr_[alternating][16] + distr_[alternating][17];
	moments[17] *= molecularVelocityCubic;	//!function for developing moment my

	moments[18] = distr_[alternating][6] - distr_[alternating][7] - distr_[alternating][8] + distr_[alternating][9] + distr_[alternating][10] - distr_[alternating][11] - distr_[alternating][12] + distr_[alternating][13];
	moments[18] *= molecularVelocityCubic;	//!function for developing moment mz
}


void Node::calcAllMoments_GSbasis(Control& bc, double* distr, double* moments, bool& alternating)
{
	double molecularVelocity = bc.getMolecularVelocity();
	double molecularVelocitySquared = pow(molecularVelocity, 2.);
	double molecularVelocityCubic = pow(molecularVelocity, 3.);
	double molecularVelocityQuartic = pow(molecularVelocity, 4.);

	//!calc moments --> save operations by hard coding
	moments[0] = 0.;
	for (int dir = 0; dir < 19; dir++) moments[0] += distr[dir];	//!density

	moments[1] = 0.;
	for (int dir = 6; dir < 18; dir++) moments[1] += distr[dir];
	moments[1] -= distr[18];
	moments[1] *= molecularVelocitySquared;	//!function for developing energy

	moments[2] = 0.;
	for (int dir = 0; dir < 6; dir++) moments[2] -= 2. * distr[dir];
	for (int dir = 6; dir < 19; dir++) moments[2] += distr[dir];
	moments[2] *= molecularVelocityQuartic;	//function for developing energy power2

	moments[3] = distr[0] - distr[2] + distr[6] + distr[7] - distr[10] - distr[11] - distr[14] + distr[15] + distr[16] - distr[17];
	moments[3] *= molecularVelocity;	//!function for developing momentum in x-direction

	moments[4] = -2. * distr[0] + 2 * distr[2] + distr[6] + distr[7] - distr[10] - distr[11] - distr[14] + distr[15] + distr[16] - distr[17];
	moments[4] *= molecularVelocityCubic;	//!function for developing heat flux in x-direction

	moments[5] = distr[1] - distr[3] + distr[8] + distr[9] - distr[12] - distr[13] + distr[14] - distr[15] + distr[16] - distr[17];
	moments[5] *= molecularVelocity;	//!function for developing momentum in y-direction

	moments[6] = -2. * distr[1] + 2 * distr[3] + distr[8] + distr[9] - distr[12] - distr[13] + distr[14] - distr[15] + distr[16] - distr[17];
	moments[6] *= molecularVelocityCubic;	//!function for developing heat flux in y-direction

	moments[7] = distr[4] - distr[5] + distr[6] - distr[7] + distr[8] - distr[9] + distr[10] - distr[11] + distr[12] - distr[13];
	moments[7] *= molecularVelocity;	//!function for developing momentum in z-direction

	moments[8] = -2. * distr[4] + 2 * distr[5] + distr[6] - distr[7] + distr[8] - distr[9] + distr[10] - distr[11] + distr[12] - distr[13];
	moments[8] *= molecularVelocityCubic;	//!function for developing heat flux in z-direction

	moments[9] = 2. * distr[0] - distr[1] + 2. * distr[2] - distr[3] - distr[4] - distr[5] + distr[6] + distr[7] - 2. * distr[8] - 2. * distr[9] + distr[10] + distr[11]
		- 2. * distr[12] - 2. * distr[13] + distr[14] + distr[15] + distr[16] + distr[17];
	moments[9] *= molecularVelocitySquared;	//!function for developing pressure xx

	moments[10] = -2. * distr[0] + distr[1] - 2. * distr[2] + distr[3] + distr[4] + distr[5] + distr[6] + distr[7] - 2. * distr[8] - 2. * distr[9] + distr[10] + distr[11]
		- 2. * distr[12] - 2. * distr[13] + distr[14] + distr[15] + distr[16] + distr[17];
	moments[10] *= molecularVelocityQuartic;	//function for developing shear stress xx

	moments[11] = distr[1] + distr[3] - distr[4] - distr[5] - distr[6] - distr[7] - distr[10] - distr[11] + distr[14] + distr[15] + distr[16] + distr[17];
	moments[11] *= molecularVelocitySquared;	//!function for developing pressure zz

	moments[12] = -distr[1] - distr[3] + distr[4] + distr[5] - distr[6] - distr[7] - distr[10] - distr[11] + distr[14] + distr[15] + distr[16] + distr[17];
	moments[12] *= molecularVelocityQuartic;	//!function for developing shear stress zz

	moments[13] = -distr[14] - distr[15] + distr[16] + distr[17];
	moments[13] *= molecularVelocitySquared; //!function for developing pressure xy

	moments[14] = distr[8] - distr[9] - distr[12] + distr[13];
	moments[14] *= molecularVelocitySquared; //!function for developing pressure yz

	moments[15] = distr[6] - distr[7] - distr[10] + distr[11];
	moments[15] *= molecularVelocitySquared; //!function for developing pressure xz

	moments[16] = -distr[6] - distr[7] + distr[10] + distr[11] - distr[14] + distr[15] + distr[16] - distr[17];
	moments[16] *= molecularVelocityCubic;	//!function for developing moment mx

	moments[17] = distr[8] + distr[9] - distr[12] - distr[13] - distr[14] + distr[15] - distr[16] + distr[17];
	moments[17] *= molecularVelocityCubic;	//!function for developing moment my

	moments[18] = distr[6] - distr[7] - distr[8] + distr[9] + distr[10] - distr[11] - distr[12] + distr[13];
	moments[18] *= molecularVelocityCubic;	//!function for developing moment mz
}

void Node::calcAllMomentsPreCol_GSbasis(Control& bc, double* moments)
{
	double molecularVelocity = bc.getMolecularVelocity();
	double molecularVelocitySquared = pow(molecularVelocity, 2.);
	double molecularVelocityCubic = pow(molecularVelocity, 3.);
	double molecularVelocityQuartic = pow(molecularVelocity, 4.);

	//!calc moments --> save operations by hard coding
	moments[0] = 0.;
	for (int dir = 0; dir < 19; dir++) moments[0] += distrPreCol_[dir];	//!density

	moments[1] = 0.;
	for (int dir = 6; dir < 18; dir++) moments[1] += distrPreCol_[dir];
	moments[1] -= distrPreCol_[18];
	moments[1] *= molecularVelocitySquared;	//!function for developing energy

	moments[2] = 0.;
	for (int dir = 0; dir < 6; dir++) moments[2] -= 2. * distrPreCol_[dir];
	for (int dir = 6; dir < 19; dir++) moments[2] += distrPreCol_[dir];
	moments[2] *= molecularVelocityQuartic;	//function for developing energy power2

	moments[3] = distrPreCol_[0] - distrPreCol_[2] + distrPreCol_[6] + distrPreCol_[7] - distrPreCol_[10] - distrPreCol_[11] - distrPreCol_[14] + distrPreCol_[15] + distrPreCol_[16] - distrPreCol_[17];
	moments[3] *= molecularVelocity;	//!function for developing momentum in x-direction

	moments[4] = -2. * distrPreCol_[0] + 2 * distrPreCol_[2] + distrPreCol_[6] + distrPreCol_[7] - distrPreCol_[10] - distrPreCol_[11] - distrPreCol_[14] + distrPreCol_[15] + distrPreCol_[16] - distrPreCol_[17];
	moments[4] *= molecularVelocityCubic;	//!function for developing heat flux in x-direction

	moments[5] = distrPreCol_[1] - distrPreCol_[3] + distrPreCol_[8] + distrPreCol_[9] - distrPreCol_[12] - distrPreCol_[13] + distrPreCol_[14] - distrPreCol_[15] + distrPreCol_[16] - distrPreCol_[17];
	moments[5] *= molecularVelocity;	//!function for developing momentum in y-direction

	moments[6] = -2. * distrPreCol_[1] + 2 * distrPreCol_[3] + distrPreCol_[8] + distrPreCol_[9] - distrPreCol_[12] - distrPreCol_[13] + distrPreCol_[14] - distrPreCol_[15] + distrPreCol_[16] - distrPreCol_[17];
	moments[6] *= molecularVelocityCubic;	//!function for developing heat flux in y-direction

	moments[7] = distrPreCol_[4] - distrPreCol_[5] + distrPreCol_[6] - distrPreCol_[7] + distrPreCol_[8] - distrPreCol_[9] + distrPreCol_[10] - distrPreCol_[11] + distrPreCol_[12] - distrPreCol_[13];
	moments[7] *= molecularVelocity;	//!function for developing momentum in z-direction

	moments[8] = -2. * distrPreCol_[4] + 2 * distrPreCol_[5] + distrPreCol_[6] - distrPreCol_[7] + distrPreCol_[8] - distrPreCol_[9] + distrPreCol_[10] - distrPreCol_[11] + distrPreCol_[12] - distrPreCol_[13];
	moments[8] *= molecularVelocityCubic;	//!function for developing heat flux in z-direction

	moments[9] = 2. * distrPreCol_[0] - distrPreCol_[1] + 2. * distrPreCol_[2] - distrPreCol_[3] - distrPreCol_[4] - distrPreCol_[5] + distrPreCol_[6] + distrPreCol_[7] - 2. * distrPreCol_[8] - 2. * distrPreCol_[9] + distrPreCol_[10] + distrPreCol_[11]
		- 2. * distrPreCol_[12] - 2. * distrPreCol_[13] + distrPreCol_[14] + distrPreCol_[15] + distrPreCol_[16] + distrPreCol_[17];
	moments[9] *= molecularVelocitySquared;	//!function for developing pressure xx

	moments[10] = -2. * distrPreCol_[0] + distrPreCol_[1] - 2. * distrPreCol_[2] + distrPreCol_[3] + distrPreCol_[4] + distrPreCol_[5] + distrPreCol_[6] + distrPreCol_[7] - 2. * distrPreCol_[8] - 2. * distrPreCol_[9] + distrPreCol_[10] + distrPreCol_[11]
		- 2. * distrPreCol_[12] - 2. * distrPreCol_[13] + distrPreCol_[14] + distrPreCol_[15] + distrPreCol_[16] + distrPreCol_[17];
	moments[10] *= molecularVelocityQuartic;	//function for developing shear stress xx

	moments[11] = distrPreCol_[1] + distrPreCol_[3] - distrPreCol_[4] - distrPreCol_[5] - distrPreCol_[6] - distrPreCol_[7] - distrPreCol_[10] - distrPreCol_[11] + distrPreCol_[14] + distrPreCol_[15] + distrPreCol_[16] + distrPreCol_[17];
	moments[11] *= molecularVelocitySquared;	//!function for developing pressure zz

	moments[12] = -distrPreCol_[1] - distrPreCol_[3] + distrPreCol_[4] + distrPreCol_[5] - distrPreCol_[6] - distrPreCol_[7] - distrPreCol_[10] - distrPreCol_[11] + distrPreCol_[14] + distrPreCol_[15] + distrPreCol_[16] + distrPreCol_[17];
	moments[12] *= molecularVelocityQuartic;	//!function for developing shear stress zz

	moments[13] = -distrPreCol_[14] - distrPreCol_[15] + distrPreCol_[16] + distrPreCol_[17];
	moments[13] *= molecularVelocitySquared; //!function for developing pressure xy

	moments[14] = distrPreCol_[8] - distrPreCol_[9] - distrPreCol_[12] + distrPreCol_[13];
	moments[14] *= molecularVelocitySquared; //!function for developing pressure yz

	moments[15] = distrPreCol_[6] - distrPreCol_[7] - distrPreCol_[10] + distrPreCol_[11];
	moments[15] *= molecularVelocitySquared; //!function for developing pressure xz

	moments[16] = -distrPreCol_[6] - distrPreCol_[7] + distrPreCol_[10] + distrPreCol_[11] - distrPreCol_[14] + distrPreCol_[15] + distrPreCol_[16] - distrPreCol_[17];
	moments[16] *= molecularVelocityCubic;	//!function for developing moment mx

	moments[17] = distrPreCol_[8] + distrPreCol_[9] - distrPreCol_[12] - distrPreCol_[13] - distrPreCol_[14] + distrPreCol_[15] - distrPreCol_[16] + distrPreCol_[17];
	moments[17] *= molecularVelocityCubic;	//!function for developing moment my

	moments[18] = distrPreCol_[6] - distrPreCol_[7] - distrPreCol_[8] + distrPreCol_[9] + distrPreCol_[10] - distrPreCol_[11] - distrPreCol_[12] + distrPreCol_[13];
	moments[18] *= molecularVelocityCubic;	//!function for developing moment mz
}
