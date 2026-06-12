// This header file is included in the main solver file airfoilLaminarFoam.C 
// to read the initial fields for pressure, velocity, and kinematic viscosity 
// from the case files.

//Including all relevant header files here:
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
#include "IOobject.H"      // IOobject class for file I/O operations
#include "IOdictionary.H"  // IOdictionary class for reading dictionaries from files

Info<< "Reading field pressure value 'P' \n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.time().name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field velocity value 'U' \n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.time().name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Reading physical properties from the physicalProperties file, which contains the kinematic viscosity value for the laminar flow case.
IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Looking up the kinematic viscosity value for the laminar flow case from the physicalProperties file, which is used in the momentum equation to model the viscous effects in the flow.
Info<< "Reading field kinematic viscosity value 'nu' \n" << endl;
dimensionedScalar nu
(
    physicalProperties.lookup("nu")
);

//Reading/Writing the face flux field 'phi', which is used in the momentum and pressure correction equations to represent the mass flux across the cell faces. The field is read from the case files and is also written to file for post-processing and visualization.
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.time().name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf() // Initializing the face flux field based on the velocity field and the mesh geometry, using the finite volume method to calculate the flux across the cell faces.
);

// Construct viscosity model from constant/physicalProperties.
// viscosityModel::New reads viscosityModel constant; nu 1.5e-05;
// and returns a constantViscosityModel registered with the mesh registry.
Info<< "Constructing viscosity model\n" << endl;
autoPtr<viscosityModel> viscosity
(
    viscosityModel::New(mesh)
);

// Construct incompressible momentum transport model.
// Passes the viscosity object as the third argument (const viscosity&).
// This registers the model in the mesh object registry under the name
// that forceCoeffs::devTau() looks up to compute viscous stress.
Info<< "Constructing momentum transport model\n" << endl;
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New
    (
        U,
        phi,
        viscosity()
    )
);