# 🌊 CFD Workspace — Einie-ai

A personal research workspace for Computational Fluid Dynamics projects, organised across three core domains.

---

## 📁 Repository Structure

```
cfd-workspace/
├── external-internal-flow/     # External & Internal Flow CFD
├── propulsion-numerics/        # Propulsion-Based Numerical Modelling
└── sciml-cfd/                  # SciML-Based CFD Numeric Case Studies
    └── project-1-pinn/         # ✅ Physics-Informed Neural Network (PINN) Model
```

---

## 🗂️ Project Categories

### 1. External & Internal Flow CFD [`/external-internal-flow`]
Classical CFD studies covering:
- External aerodynamics (bluff bodies, aerofoils, boundary layers)
- Internal flows (pipe flow, channel flow, duct networks)
- Turbulence modelling (RANS, LES, DNS)

### 2. Propulsion-Based Numerical Modelling [`/propulsion-numerics`]
Numerical simulation focused on propulsion systems:
- Nozzle flow & shock structures
- Combustion and reacting flows
- Rocket / jet engine component modelling

### 3. SciML CFD Case Studies [`/sciml-cfd`]
Scientific Machine Learning applied to fluid dynamics:
- Physics-Informed Neural Networks (PINNs)
- Neural operator methods (FNO, DeepONet)
- Data-driven turbulence closure

---

## ✅ Active Projects

| # | Name | Category | Status |
|---|------|----------|--------|
| 1 | [PINN Baseline Model](./sciml-cfd/project-1-pinn/) | SciML CFD | 🟢 Active |

---

## ⚙️ Environment & Tooling

| Tool | Purpose |
|------|---------|
| Python 3.10+ | Primary language |
| PyTorch / JAX | ML frameworks |
| OpenFOAM | Classical CFD solver |
| NumPy / SciPy | Numerical computing |
| Matplotlib / ParaView | Visualisation |

---

## 📜 Conventions

- Each project lives in its own sub-directory with its own `README.md`
- Results and large datasets go in a `results/` folder (`.gitignore`d)
- Use `requirements.txt` or `environment.yml` per project

---

*Maintained by [@Einie-ai](https://github.com/Einie-ai)*
