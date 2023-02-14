/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    fsinetFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving netting.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   This new solver can solve moving nets with various solidities and shapes. The nets is
   modelled by a dynamic thin porous media zone where the position of the porous zone is defined
   by a external file and the values of source term is calculated based on the forces on the respresented nets.
   
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

#include "netPanel.H"
#include "OFstream.H"
#include "Pstream.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    OFstream* myOutFilePtr = NULL;

    if (Pstream::master())
    {
        // Open the file
        myOutFilePtr = new OFstream("velocity_on_elements.txt");
    }
    
while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        // read from outside
        // need to confirm ... might be wrong data structure
        Info<< ">>> Start fsi function..."<< "  ClockTime = " << runTime.elapsedClockTime() << " s"<<endl;
        while (not exists("./constant/position_flag")) {
            sleep(0.1);
            Info<<">>>    Waiting 0.01s for posi"<<endl;
        }
        
        Foam::IOdictionary structuralPositions(
                IOobject(
                        "posi",
                        runTime.constant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE));
        
        Nettings.readPosi(structuralPositions);
        
        while (not exists("./constant/fh_flag")){
            sleep(0.1);
            Info<<">>>    Waiting 0.01s for fh"<<endl;
        }

        IOdictionary structuralFh(
                IOobject(
                        "Fh",
                        runTime.constant(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE));
        
        Nettings.readForce(runTime.value(),structuralFh);

        List<pointField> gatheredU(numberP);
        gatheredU[Pstream::myProcNo()] = pointField(U);
        Pstream::gatherList(gatheredU);

        
        Nettings.updateVelocity(gatheredU,gatheredCentres,thresholdLength);
        

//        Info<<"velocity on the center of net panels are \n"<<Nettings.FluidU()<<endl;
        // write the Nettings.fluidVelocity(); to a extrinal files
//        Info<< "Start writing velocity"<< "  ClockTime = " << runTime.elapsedClockTime() << " s"<<endl;
        if (Pstream::master())
        {
            OFstream& myOutFile = *myOutFilePtr;
            myOutFile
                    << "The velocities at " << runTime.timeName()<< " s are: "<< Nettings.FluidU()  << endl;
        }
//        Info<< "Finish writing velocity"<< "  ClockTime = " << runTime.elapsedClockTime() << " s"<<endl;

        Info<< ">>> Finish fsi function..."<< "  ClockTime = " << runTime.elapsedClockTime() << " s\n "<<endl;
        runTime.write();
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }


    Info<< "End\n" << endl;


    return 0;
}

// ************************************************************************* //