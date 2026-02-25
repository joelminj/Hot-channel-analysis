# PWR Hot Channel Thermal-Hydraulics Analysis (MATLAB)

## Overview

This project implements a one-dimensional, steady-state thermal-hydraulics model for the **hot channel of a Pressurized Water Reactor (PWR)** using MATLAB.

The objective is to evaluate:

- Axial coolant enthalpy rise
- Convective heat transfer at the cladding surface
- Radial temperature distribution in the fuel rod
- Onset of nucleate boiling conditions

The model resolves the temperature profile from the coolant bulk to the fuel centerline under nominal operating conditions.

---

## Physical Model

The active core height is discretized into **N axial nodes**, and conservation of energy is applied at each node.

The following physical phenomena are modeled:

### Axial Power Distribution

The axial heat generation rate is calculated using a chopped cosine power profile:

$$
q_i = \bar{q}'_{hot} L_{node} F_{z,i}
$$

Surface heat flux:

$$
q''_i =
\frac{q_i}{2 \pi r_{co} L_{node}}
$$

---

### Coolant Energy Balance

The coolant enthalpy rise is calculated using:

$$
h_{out,i} =
h_{in,i}
+
\frac{q_i}{\dot{m}_{ch}}
$$

Fluid properties are evaluated at bulk coolant temperature.

---

### Convective Heat Transfer

Two correlations are implemented:

#### Dittus-Boelter correlation

$$
Nu =
0.023 Re^{0.8} Pr^{0.4}
$$

$$
h =
\frac{Nu k}{D_h}
$$

#### Bernath correlation

Applicable for subcooled and extended boiling regimes:

$$
h =
10890
\left(
\frac{D_e}{D_e + D_i}
\right)
+
\frac{48 V}{D_e^{0.6}}
$$

The limiting heat transfer coefficient is used.

---

### Radial Heat Conduction

Radial conduction is solved across:

- Cladding
- Helium gap
- Fuel pellet

Cladding inner temperature:

$$
T_{ci} =
T_{co}
+
\frac{q_i \ln(r_{co}/r_{ci})}
{2 \pi k_{clad} L}
$$

Fuel outer temperature:

$$
T_{fo} =
T_{ci}
+
\frac{q_i \ln(r_{ci}/r_{fo})}
{2 \pi k_{gap} L}
$$

Fuel centerline temperature:

$$
T_c =
T_{fo}
+
\frac{q_i}
{4 \pi k_{fuel} L}
$$

---

### Nucleate Boiling Limit

Thom correlation is used:

$$
T_{thom}
=
T_{sat}
+
22.5
\sqrt{q''}
\exp
\left(
-\frac{P}{8.7}
\right)
$$

This allows prediction of subcooled boiling onset.

---

## Features

- Axial coolant temperature distribution
- Fuel centerline temperature calculation
- Cladding temperature prediction
- Heat transfer coefficient evaluation
- Boiling onset prediction
- Multiple heat transfer correlations

---

## Code Structure

hot_rod.m


Main script performing thermal-hydraulic calculations.

---

## Requirements

Software:

- MATLAB R2020 or newer

Optional:

- XSteam or steam table library

---

## How to Run

Run the script:

Open the script in matlab and run
