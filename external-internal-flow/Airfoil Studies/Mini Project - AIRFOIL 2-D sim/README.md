# External & Internal Flow CFD

Mini-Project: Built a standalone, custom-compiled CFD solver from scratch within the OpenFOAM v13 framework for simulating 2D incompressible laminar flow over a NACA 0012 airfoil — without inheriting from any existing solver module.
The solver implements the full four-step SIMPLE algorithm using raw fvc::/fvm:: operators: momentum predictor, HbyA extraction, pressure Poisson solve via GAMG + GaussSeidel, and explicit velocity correction. A bespoke createFields.H was authored to handle field initialisation, viscosity registration, and interfacing with the forceCoeffs functionObject for aerodynamic coefficient extraction (Cl, Cd, Cm).
The test case ran at Re ≈ 2×10⁶ (U∞ = 30 m/s, ν = 1.5×10⁻⁵ m²/s) on the airFoil2D polyMesh (~12,000 hexahedral cells), completing 10 SIMPLE iterations in 0.75 s CPU time. While the velocity field produced qualitatively correct flow topology — leading edge acceleration, no-slip wall layer, and wake deficit — solver divergence was observed in the force coefficients by iteration 10 (Cd = −2.75), attributable to a geometry mismatch in the baseline mesh, incomplete Rhie-Chow interpolation, overly conservative relaxation factors, and the use of a laminar model at a physically turbulent Reynolds number.
The exercise demonstrated end-to-end ownership of a custom OpenFOAM application: C++ solver design, wmake build configuration, case setup, and quantitative post-processing in ParaView — with a rigorous root-cause analysis of the numerical limitations introduced by deliberate solver simplifications.

**Still in progress to improve convergence over larger number of time steps.

Tools: OpenFOAM v13 · C++ · ParaView · wmake · SIMPLE Algorithm

