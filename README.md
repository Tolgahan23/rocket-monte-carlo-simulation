# Stochastic Flight Dynamics: Rocket Trajectory Monte Carlo Simulation

This project is a MATLAB-based simulation designed to analyze the dispersion and reliability of a rocket's trajectory under stochastic atmospheric conditions [1]. It utilizes the Monte Carlo method to execute 2,000 independent iterations, accounting for real-world thermodynamic variability [1, 2].

## 🚀 Project Overview
In aerospace engineering, deterministic models often fail to capture the uncertainty of the environment, assuming a "perfect" and static atmosphere [3]. This simulation introduces randomness into key atmospheric parameters to move beyond a single point estimate and create a comprehensive probabilistic landscape [4, 5]. By doing so, it visualizes the statistically supported "Impact Dispersion Zone" and quantifies the probability of mission success [4].

## ✨ Key Features
* **Stochastic Modeling:** Incorporates random variables for initial sea-level temperature (T0) and pressure (P0) using Normal (Gaussian) distributions to simulate natural tropospheric noise [2, 5].
* **Kinematic & Atmospheric Analysis:** Solves 2D equations of motion for a variable-mass system with a precise discrete time-step of dt = 0.01s [6, 7]. It seamlessly integrates the International Standard Atmosphere (ISA) equations to dynamically calculate changing air density and drag force [8, 9].
* **Automated Data Visualization:** Automatically aggregates data from all 2,000 flight cycles to generate Probability Density Functions (PDF), Cumulative Distribution Functions (CDF) for defining safety zones, and 2D Density Maps (heatmaps) to identify the impact of environmental uncertainty [5, 10, 11].

## 🛠 Tech Stack
* **Language:** MATLAB [10].
* **Methodology:** Monte Carlo Simulation, Discrete Time-Step Numerical Integration, Probabilistic Modeling & Risk Assessment [4-6].

## 📊 How to Run
1. Clone the repository: `git clone https://github.com/Tolgahan23/rocket-monte-carlo-simulation.git`
2. Open the project folder in MATLAB [10].
3. Run the main script to initiate the 2,000 flight iterations. The tool will automatically perform post-processing and display the trajectory dispersion plots, PDFs, CDFs, and 2D Maps without requiring manual data analysis [10].

## 📊 Simulation Results
<img width="1060" height="860" alt="Figure_1" src="https://github.com/user-attachments/assets/254e8e1f-2f31-4faa-8299-bfe27ecade8e" />
<img width="1064" height="861" alt="Figure_2" src="https://github.com/user-attachments/assets/de121549-c3c9-4485-9fe9-1251f1b801c3" />
<img width="1180" height="849" alt="Figure_9" src="https://github.com/user-attachments/assets/f7a396ce-1ad7-4ab7-9521-88d613164002" />
