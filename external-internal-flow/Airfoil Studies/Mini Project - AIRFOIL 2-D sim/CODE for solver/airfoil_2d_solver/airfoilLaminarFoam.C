/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
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
    airfoilLaminarFoam

Description
Solves for the Velocity, pressure and the phi values for a laminar flow over a 2-D airfoil. 
Is slightly/amateurely based on te simple algorithm.

\*---------------------------------------------------------------------------*/

#include "argList.H"
// fvCFD.H contents below:
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"         // runTime, time loop control
#include "fvMesh.H"       // mesh geometry + cell/face access
#include "fvc.H"          // explicit operators: fvc::div, fvc::grad, fvc::laplacian
#include "fvMatrices.H"   // fvMatrix<Type>: the discrete equation matrix
#include "fvm.H"          // implicit operators: fvm::ddt, fvm::div, fvm::laplacian
#include "linear.H"       // linear interpolation
#include "uniformDimensionedFields.H"
#include "calculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "adjustPhi.H"    // flux correction at boundaries
#include "findRefCell.H"  // pressure reference cell utility
#include "fvModels.H"       // turbulence model base class
#include "fvConstraints.H"   // constraint utilities for Rhie-Chow flux correction and pressure boundary enforcement
#include "airfoilLaminarFoam.H" // Custom header file for the standalone laminar SIMPLE airfoil solver, which pulls in only the extra special classes actually used in the solver.

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    

    // Creating the simpleControl object which owns the residual control dictionary in which residuals of all fields are documented.
    simpleControl simple(mesh);

    Info<<"\n Reading the initial fields... \n"<<nl<<endl;
    #include "createFields.H"

    // Printing the basic initial field values for the case, which are read from the createFields.H file.
    Info<<"\n Printing basic initial field values for the case: \n"<<nl<<endl;
    Info<<"Initial pressure field: p = "<<gAverage(p.internalField())<<nl<<endl;
    Info<<"Initial velocity field: U = "<<gAverage(U.internalField())<<nl<<endl;
    Info<<"Initial kinematic viscosity: nu = "<<nu<<nl<<endl;

    // Looking up the reference cell and value for pressure from the fvSolution file, which is used to set a reference point for the pressure field to ensure a unique solution.
    label pRefCell = 0; // Default reference cell index for pressure, which can be overridden by the value specified in the fvSolution file. The reference cell is used to set a reference point for the pressure field to ensure a unique solution, as pressure is defined up to an arbitrary constant.
    scalar pRefValue = 0.0; // Default reference value for pressure, which can be overridden by the value specified in the fvSolution file. The reference value is used to set a reference point for the pressure field to ensure a unique solution, as pressure is defined up to an arbitrary constant.
    setRefCell(p,simple.dict(),pRefCell, pRefValue); // Setting the reference cell and value for pressure based on the values specified in the simpleControl dictionary, which is read from the fvSolution file. The setRefCell function is used to set the reference cell and value for the pressure field, ensuring that the pressure field has a unique solution.
    

    Info<<"\n Reading the reference cell and value for pressure from the fvSolution file... \n"<<nl<<endl;
    Info<<"Reference cell for pressure: pRefCell = "<<pRefCell<<nl<<endl;
    Info<<"Reference value for pressure: pRefValue = "<<pRefValue<<nl<<endl;

    Info<<"\n Starting the SIMPLE based time loop... \n"<<nl<<endl;

    while(simple.loop(runTime))
    {
        Info<<" Time Iteration: "<<runTime.time().value()<<nl<<endl;

        // SIMPLE based alogorithm starts:

        // 1. Momentum equation is discretized:

        fvVectorMatrix Ueqn
        (
            fvm::div(phi,U) - fvm::laplacian(nu,U)
        );
        Ueqn.relax(); // Relaxing the momentum equation to enhance stability of the solution. 
        // The relaxation factor is defined in the fvSolution file. 
        // Here the diagonal elements of the [A] are inflated for better conevergence of the solution.
        
        if (simple.momentumPredictor()) // checking if the momentum predictor step is enabled in the simpleControl dictionary. 
        // If it is enabled, the momentum equation is solved before the pressure correction step. 
        // This can help to improve convergence in some cases, especially when the flow is highly unsteady 
        // or when there are strong velocity gradients.
        {
            solve(Ueqn == -fvc::grad(p)); // Solving the momentum equation for the velocity field. 
        }
    
        //Extracting [A] and [H] vector values from the momentum equation:
        volScalarField A= Ueqn.A(); // Diagonal coefficients of the matrix [A] in the momentum equation.
        volVectorField H = Ueqn.H(); // Off-diagonal coefficients of the matrix [A] in the momentum equation, which are stored in the vector field H.

        //2. Pressure Correction equation is discretized:
        volScalarField A_inv(1.0/A); // Inverse of the diagonal coefficients of the matrix [A] in the momentum equation, used in the pressure correction step.
        volVectorField HbyA(constrainHbyA(A_inv * H, U, p)); // The product of the inverse of the diagonal coefficients and the off-diagonal coefficients, used in the pressure correction step.
        surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) & mesh.Sf()); // The flux of the HbyA field across the cell faces, used in the pressure correction step (divergence function)
        adjustPhi(phiHbyA, U, p); // Adjusting the flux field phibyA to ensure mass conservation at the boundaries, based on the current velocity and pressure fields at the cell faces and centres.
        
        //Pressure equation is solved for the pressure correction field, which is used to correct the velocity field to ensure mass conservation.
        fvScalarMatrix pEqn
        (
             fvm::laplacian(A_inv,p)== fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue); // Setting a reference cell for the pressure field to ensure a unique solution, as pressure is defined up to an arbitrary constant.
        pEqn.solve(); // Solving the pressure correction equation for the pressure field.
        
        // 4. Flux field is corrected based on the corrected velocity field:
        phi= phiHbyA - pEqn.flux(); // Correcting the flux field based on the corrected velocity field, ensuring that the flux field is consistent with the corrected velocity field and the mesh geometry.
    
        
        p.relax(); // Relaxing the pressure field to enhance stability of the solution. The relaxation factor is defined in the fvSolution file. Here the pressure field is relaxed to improve convergence of the solution.
        constrainPressure(p, U, phiHbyA, A_inv); // Constraining the pressure field to ensure that the pressure correction is consistent with the velocity and flux fields at the boundary face conditions, based on the current velocity and flux fields and the coefficients from the momentum equation.

        // 3. Velocity field is corrected based on the pressure correction:
        U= HbyA - A_inv*fvc::grad(p); // Correcting the velocity field based on the pressure correction, ensuring that the corrected velocity field satisfies the momentum equation and mass conservation.
        U.correctBoundaryConditions(); // Correcting the velocity field at the boundaries to ensure that the corrected velocity field satisfies the boundary conditions, based on the current velocity and pressure fields and the coefficients from the momentum equation.

        
       // Writing the corrected fields to file for post-processing and visualization:
        runTime.write(); // Writing the current time step to file, which is used for time control and post-processing.
        
        Info<<"\n"<<nl<<endl;
    }
    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
