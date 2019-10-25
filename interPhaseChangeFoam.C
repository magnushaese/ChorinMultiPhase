/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    ReboudInterPhaseChangeFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids with phase-change
    (e.g. cavitation).  Uses a VOF (volume of fluid) phase-fraction based
    interface capturing approach. Includes a correction for the 
    turbulent viscosity in the two phase region due to Reboud.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate cavitation
    but other mechanisms of phase-change are supported within this solver
    framework.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pimpleControl pimple(mesh);

#   include "readGravitationalAcceleration.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "createTimeControls.H"
#   include "correctPhi.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"
#   include "readRebExp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

// solve for volume fraction
        #include "alphaEqnSubCycle.H"
        twoPhaseProperties->correct();


//// FRACTIONAL STEP	
// STEP 1 -- velocity prediction
//
	Uold = U;
        phi = (fvc::interpolate(U) & mesh.Sf());
        dU = dt*(fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi,U)); //rho^n+1 already inside fvc::[...]
	U = dU + rho.oldTime()/rho*Uold;
//	#include "continuityErrs.H"
	U.correctBoundaryConditions();
//	phi = (fvc::interpolate(U) & mesh.Sf());

// STEP 2 -- implicit pressure calculation from Poisson equation
//
       Pair<tmp<volScalarField> > vDotP = twoPhaseProperties->vDotP(); //get cavitation source terms for present (n+1) velocity divergence
       const volScalarField& vDotcP = vDotP[0]();
       const volScalarField& vDotvP = vDotP[1]();

	fvScalarMatrix pdEqn
        (
                fvm::laplacian(pd)
//              + (vDotvP - vDotcP)*(rho*gh - pSat) + fvm::Sp(vDotvP - vDotcP, pd)
	);
        pdEqn.setReference(pdRefCell, pdRefValue);

	solve(pdEqn == (scalar(1.0)/dt*rho)*(fvc::div(U) + (vDotvP - vDotcP)*(rho*gh - pSat) + fvc::Sp(vDotvP - vDotcP, pd)));

// STEP 3 -- final velocity correction step
	U -= dt/rho*fvc::grad(pd); 
	U.correctBoundaryConditions();
//	phi = (fvc::interpolate(U) & mesh.Sf());
// END FRACTIONAL STEP

        turbulence->correct();

	runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
