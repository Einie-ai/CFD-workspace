/*----------------------------------------------------------*\
  airfoilLaminarFoam.H
  Custom include header for the standalone laminar SIMPLE
  airfoil solver.  Pulls in only the classes actually used:
    - fvCFD.H              : mesh, time, fvc/fvm operators,
                             fvMatrix, field types
    - simpleControl.H      : SIMPLE loop convergence object
    - singlePhaseTransportModel.H : constant-nu viscosity
    - adjustPhi.H          : phi boundary correction
    - constrainHbyA.H      : Rhie-Chow velocity constraint
    - constrainPressure.H  : pressure boundary enforcement
    - findRefCell.H        : reference cell/value for p
\*----------------------------------------------------------*/

#ifndef airfoilLaminarFoam_H
#define airfoilLaminarFoam_H

#include "simpleControl.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"

// OF13 momentum transport model — replaces singlePhaseTransportModel
// which does not exist in OF13. Reads nu from physicalProperties directly.
#include "incompressibleMomentumTransportModel.H"

#include "viscosity.H" // For the nu() function, which is used to read the kinematic viscosity value from the physicalProperties file.
#include "viscosityModel.H"

#endif