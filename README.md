```markdown
# ⚡ Effective Mass Calculation using Quantum Espresso (thermo_pw.x)

This repository provides a **simple workflow** to calculate the **effective mass of electrons and holes** using:

- Quantum Espresso (QE)
- thermo_pw.x (post-processing)
- Python (optional for plotting)

---

## 🎯 Overview

Effective mass describes how electrons and holes respond to external forces and is essential for:

- Carrier transport  
- Mobility  
- Semiconductor properties  
- Electronic device performance  

---

## ⚙️ Workflow

```

SCF Calculation
↓
NSCF Calculation
↓
thermo\_pw\.x Analysis
↓
Effective Mass Extraction

````

---

## 📁 Files Included

- `scf.in` → SCF input file  
- `nscf.in` → NSCF input file   
- `thermo_pw.in` → thermo_pw.x input file  
- `bands.dat` → Band data  
- `effective_mass.py` → Optional Python script  

---

## 🔧 Step-by-Step

### 1️⃣ SCF Calculation

```bash
pw.x < scf.in > scf.out
````

***

### 2️⃣ NSCF Calculation

```bash
pw.x < nscf.in > nscf.out
```

***

### 3️⃣ thermo\_pw\.x Analysis

```bash
thermo_pw.x < thermo_pw.in > thermo_pw.out
```

👉 thermo\_pw\.x can:

* Analyze band structure
* Extract curvature
* Assist in effective mass analysis

***

## 🧠 Effective Mass Formula

$$
m^* = \hbar^2 \left(\frac{d^2E}{dk^2}\right)^{-1}
$$

👉 Derived from **band curvature near CBM/VBM**

***

## 📊 Key Concepts

* CBM → electron effective mass
* VBM → hole effective mass
* Flat band → high mass
* Curved band → low mass

***

## 🐍 Optional Python (Curve Fitting)

```python
import numpy as np

k = np.array([...])
E = np.array([...])

coeff = np.polyfit(k, E, 2)
d2E = 2 * coeff[0]

m_eff = 1 / d2E

print("Effective Mass:", m_eff)
```

***

## 🚀 Tips

✅ Use dense k-points near CBM/VBM  
✅ Fit only small k-range  
✅ Check band structure carefully  
✅ Use consistent units

***

## 🎥 YouTube Tutorial

👉 <https://www.youtube.com/@DeobratQMatX>

***

## 📌 Requirements

* Quantum Espresso (QE)
* thermo\_pw\.x
* Python (optional)

***

## 📬 Contact

📧 <deobratqmatx@gmail.com>

***
