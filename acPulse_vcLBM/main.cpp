#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>
#include "Control.h"
#include "Lattice.h"
#include "Node.h"

using namespace std;
using namespace chrono;

ostream&
display(ostream& os, nanoseconds ns)
{
    typedef duration<int, ratio<86400>> days;
    char fill = os.fill();
    os.fill('0');
    auto d = duration_cast<days>(ns);
    ns -= d;
    auto h = duration_cast<hours>(ns);
    ns -= h;
    auto m = duration_cast<minutes>(ns);
    ns -= m;
    auto s = duration_cast<seconds>(ns);
    os << setw(2) << d.count() << "d:"
        << setw(2) << h.count() << "h:"
        << setw(2) << m.count() << "m:"
        << setw(2) << s.count() << 's' << endl << endl;
    os.fill(fill);
    return os;
};

int main(int argc, char* argv[])
{
    steady_clock::time_point beginTotal = steady_clock::now();
    steady_clock::time_point endTotal;
    steady_clock::time_point beginCalc;
    steady_clock::time_point endCalc;
    nanoseconds calcTime;
    string inputFile = "/absolute-path-to-BC-file/BoundaryConditions.txt"; //!name of input file; could also be replaced by argv for dynamic file names
    string initTimeStepString;
    Control boundaryConditions(inputFile); //!create instance of control class and read parameters from file
    Lattice simulationCase(boundaryConditions); //!create instance of lattice
    int time = 0;
    int initTime = 0;

    int writeInterval = boundaryConditions.getWriteInterval(); //!get write interval
    int maxTimeStep = boundaryConditions.getTimeStepMax();
    int maxLevel = boundaryConditions.getLevelMax();
    boundaryConditions.setTime_(static_cast<double>(time));
    bool alternatingL0 = true; //!variable for switching between distribution arrays
    bool negAlternatingL0 = !alternatingL0;
    bool alternatingL1 = true;

    time = 0;

    //!write to console
    cout << "LB acoustic pulse solver" << endl;

    cout << "initializing from uniform zero-velocity equilibrium distribution" << endl << endl;

    if (maxLevel == 0)
    {
        cout << "all parameters are being displayed in SI units" << endl << endl;
        cout << "result file(s): " << boundaryConditions.getCaseName() << "_t.vtk" << endl << endl;
        cout << "        dimensionless collision frequency = " << boundaryConditions.getCollisionFrequency() << endl;
        cout << "                           time step size = " << boundaryConditions.getTimeStep() << endl;
        cout << "                                  spacing = " << boundaryConditions.getSpacing() << endl;
        cout << "                          REYNOLDS number = " << boundaryConditions.getReynoldsNumber() << endl;
        cout << "                      kinematic viscosity = " << boundaryConditions.getKinematicViscosity() << endl;
        cout << "                                  density = " << boundaryConditions.getDensity() << endl;
        cout << "                           speed of sound = " << boundaryConditions.getSpeedOfSound() << endl;
        cout << "                              MACH number = " << boundaryConditions.getMachNumber() << endl;
        cout << "                     order of equilibrium = " << boundaryConditions.getOrderOfEquilibrium() << endl;
        cout << "                          collision model = " << boundaryConditions.getCollisionModel() << endl;
        if(boundaryConditions.getCollisionModel() == "HRR")
        {
            cout << "                 hybridization parameter = " << boundaryConditions.getHybridParameter() << endl;
        }
        cout << "       non-reflecting boundary conditions = " << boundaryConditions.getReflectionHandling() << endl;
        if(boundaryConditions.getReflectionHandling() == "yes")
        {
            cout << "                        transverse terms = " << boundaryConditions.getTransverseInclusion() << endl;
            cout << "                    K1 relaxation factor = " << boundaryConditions.getLODIK1relaxation() << endl;
        }
        cout << "                         cubic correction = " << boundaryConditions.getCubicMachCorrection() << endl;
        cout << "                            domain size x = " << boundaryConditions.getChannelSizeX_() << endl;
        cout << "                            domain size y = " << boundaryConditions.getChannelSizeY_() << endl;
        cout << "                            domain size z = " << boundaryConditions.getChannelSizeZ_() << endl;
        cout << "                     number of time steps = " << boundaryConditions.getTimeStepMax() << endl;
        cout << "                          write intervall = " << boundaryConditions.getWriteInterval() << endl << endl;

        //!initialize nodes
        simulationCase.initialization(boundaryConditions);

        //!set calculation clock
        beginCalc = steady_clock::now();

        //!print current time to console
        auto startTime = system_clock::now();
        time_t cur_time = system_clock::to_time_t(startTime);

        cout << "current date and time is: " << ctime(&cur_time) << endl;

        while (time < maxTimeStep) //!loop over all time steps
        {
            if (time == 0)
            {
                //!write initialized fields
                string file = boundaryConditions.getCaseName(); //!get case name
                ostringstream timeStep; //!create string stream for conversion from integer to string
                timeStep << time; //!convert time (integer) to string
                file.append("_" + timeStep.str()); //!add current time step to file name
                simulationCase.writeResultsVTK(file, alternatingL0, alternatingL1); //!write result file
            }

            simulationCase.solve(alternatingL0); //!solve case for this time step

            time++; //!increment time
            boundaryConditions.setTime_(static_cast<double>(time));

            if (time - initTime == 1)
            {
                endCalc = steady_clock::now();
                calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * (maxTimeStep - initTime));

                cout << "estimated total computation time: ";
                display(cout, calcTime);
            } //!end if
            else if ((time % writeInterval) == 0) //!write results every write interval
            {
                //!write result file
                string file = boundaryConditions.getCaseName(); //!get case name
                ostringstream timeStep; //!create string stream for conversion from integer to string
                timeStep << time; //!convert time (integer) to string
                file.append("_" + timeStep.str()); //!add current time step to file name
                simulationCase.writeResultsVTK(file, negAlternatingL0, alternatingL1); //!write result file

                //!write status to console
                cout << "solved " << fixed << setprecision(2) << (static_cast<double>(time) - static_cast<double>(initTime)) / (maxTimeStep - static_cast<double>(initTime)) * 100 << "% of the case" << endl;

                if (time - initTime == writeInterval)
                {
                    endCalc = steady_clock::now();
                    calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * ((maxTimeStep - initTime) / writeInterval));

                    cout << endl << "corrected computation time estimation (including time for file creation): ";
                    display(cout, calcTime);
                } //!end if
            } //!end else if
            alternatingL0 = !alternatingL0; //!switch distribution array
            negAlternatingL0 = !negAlternatingL0;
        } //!end time
    } //!end if
    else
    {
        cout << "all parameters are being displayed in SI units" << endl << endl;
        cout << "result file(s): " << boundaryConditions.getCaseName() << "_t.vtk" << endl << endl;
        cout << "dimensionless collision frequency on coarse grid = " << boundaryConditions.getCollisionFrequency() << endl;
        cout << "  dimensionless collision frequency on fine grid = " << 1. / ((2. / boundaryConditions.getCollisionFrequency()) - 0.5) << endl;
        cout << "                   time step size on coarse grid = " << boundaryConditions.getTimeStep() << endl;
        cout << "                     time step size on fine grid = " << boundaryConditions.getTimeStep() / 2. << endl;
        cout << "                          spacing on coarse grid = " << boundaryConditions.getSpacing() << endl;
        cout << "                            spacing on fine grid = " << boundaryConditions.getSpacing() / 2. << endl;
        cout << "                                 REYNOLDS number = " << boundaryConditions.getReynoldsNumber() << endl;
        cout << "                             kinematic viscosity = " << boundaryConditions.getKinematicViscosity() << endl;
        cout << "                                         density = " << boundaryConditions.getDensity() << endl;
        cout << "                                  speed of sound = " << boundaryConditions.getSpeedOfSound() << endl;
        cout << "                                     MACH number = " << boundaryConditions.getMachNumber() << endl;
        cout << "                            order of equilibrium = " << boundaryConditions.getOrderOfEquilibrium() << endl;
        cout << "                                 collision model = " << boundaryConditions.getCollisionModel() << endl;
        if(boundaryConditions.getCollisionModel() == "HRR")
        {
            cout << "                         hybridization parameter = " << boundaryConditions.getHybridParameter() << endl;
        }
        cout << "                                cubic correction = " << boundaryConditions.getCubicMachCorrection() << endl;
        cout << "              non-reflecting boundary conditions = " << boundaryConditions.getReflectionHandling() << endl;
        if(boundaryConditions.getReflectionHandling() == "yes")
        {
            cout << "                                transverse terms = " << boundaryConditions.getTransverseInclusion() << endl;
            cout << "                            K1 relaxation factor = " << boundaryConditions.getLODIK1relaxation() << endl;
        }
        cout << "                               coupling approach = " << boundaryConditions.getCouplingProcedure() << endl;
        cout << "          fine links to be replaced at interface = " << boundaryConditions.getCouplingReplacementC2F() << endl;
        cout << "        coarse links to be replaced at interface = " << boundaryConditions.getCouplingReplacementF2C() << endl;
        cout << "                     hanging node reconstruction = " << boundaryConditions.getHangingNodeReconstruction() << endl;
        cout << "                                   f2c filtering = " << boundaryConditions.getF2Cfilter() << endl;
        cout << "                                   scaling style = " << boundaryConditions.getScalingStyle() << endl;
        cout << "                                   domain size x = " << boundaryConditions.getChannelSizeX_() << endl;
        cout << "                                   domain size y = " << boundaryConditions.getChannelSizeY_() << endl;
        cout << "                                   domain size z = " << boundaryConditions.getChannelSizeZ_() << endl;
        cout << "                            number of time steps = " << boundaryConditions.getTimeStepMax() << endl;
        cout << "                                 write intervall = " << boundaryConditions.getWriteInterval() << endl;
        cout << "refinement layers at Y and Z boundaries and back = " << boundaryConditions.getRefinementLayer() << endl;
        cout << "                                         overlap = " << boundaryConditions.getOverlap() << endl;
        cout << "                    refinement layers at X front = " << boundaryConditions.getRefinementLayerX() << endl << endl;

        //initialize nodes
        simulationCase.initialization(boundaryConditions);

        //!set calculation clock
        beginCalc = steady_clock::now();

        //!print current time to console
        auto startTime = std::chrono::system_clock::now();
        time_t cur_time = system_clock::to_time_t(startTime);

        cout << "current date and time is: ";
        cout << ctime(&cur_time) << endl;

        while (time < maxTimeStep) //!loop over all time steps
        {
            if (time == 0)
            {
                //!write initialized fields
                string file = boundaryConditions.getCaseName(); //!get case name
                ostringstream timeStep; //!create string stream for conversion from integer to string
                timeStep << time; //!convert time (integer) to string
                file.append("_" + timeStep.str()); //!add current time step to file name
                simulationCase.writeResultsVTK(file, alternatingL0, alternatingL1); //!write result file
            }

            simulationCase.solve(alternatingL0, alternatingL1); //!solve case for this time step

            time++; //!increment time
            boundaryConditions.setTime_(static_cast<double>(time));

            if (time - initTime == 1)
            {
                endCalc = steady_clock::now();
                calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * (maxTimeStep - initTime));

                cout << "estimated total computation time: ";
                display(cout, calcTime);
            } //!end if
            else if ((time % writeInterval) == 0) //!write results every write interval
            {
                //!write result file
                string file = boundaryConditions.getCaseName(); //!get case name
                ostringstream timeStep; //!create string stream for conversion from integer to string
                timeStep << time; //!convert time (integer) to string
                file.append("_" + timeStep.str()); //!add current time step to file name
                simulationCase.writeResultsVTK(file, negAlternatingL0, alternatingL1); //!write result file

                //!write status to console
                cout << "solved " << fixed << setprecision(2) << (static_cast<double>(time) - static_cast<double>(initTime)) / (maxTimeStep - static_cast<double>(initTime)) * 100 << "% of the case" << endl;

                if (time - initTime == writeInterval)
                {
                    endCalc = steady_clock::now();
                    calcTime = duration_cast<nanoseconds>((endCalc - beginCalc) * ((maxTimeStep - initTime) / writeInterval));

                    cout << endl << "corrected computation time estimation (including time for file creation): ";
                    display(cout, calcTime);
                } //!end if
            } //!end else if

            alternatingL0 = !alternatingL0; //!switch distribution array
            negAlternatingL0 = !negAlternatingL0;
        } //!end time
    } //!end else

    endTotal = steady_clock::now();
    cout << endl << "given maximum number of time steps reached" << endl << endl;
    cout << "total computation time: ";
    display(cout, duration_cast<nanoseconds>(endTotal - beginTotal));
} //!end main
