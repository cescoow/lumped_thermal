# MATLAB Thermal Toolbox - GUI Documentation

### April 2024

---

## Introduction
Created by Francesco Lena in October 2024, this repository hosts MATLAB tools for simulating and analyzing heat transfer in interconnected systems. The toolbox includes two main applications:

1. **View Factor Calculator**: Uses a Monte Carlo ray-tracing method to estimate radiative view factors between surfaces loaded from STL files. This is ideal for studying systems involving radiative heat exchange.
2. **Lumped Elements Thermal Simulation**: Models temperature dynamics in systems of lumped thermal elements, accounting for conduction, radiation, and external heat sources.

These tools are experimental, designed primarily for academic and research use where customization is essential. Users should validate results against other methodologies for comprehensive analysis. For more information, see the project repository at [GitHub - Lumped Thermal](https://github.com/cescoow/lumped_thermal).

---

## View Factor Calculator Application - Comprehensive Guide

### Overview
The **View Factor Calculator** employs a Monte Carlo method to compute view factors between 3D surfaces, ideal for radiative heat transfer analysis. It uses ray tracing between an emissive surface (Body A) and a target surface (Body B), with optional obstructions. This application is meant for scientific calculations in thermal radiation and computational heat transfer studies. Validation with established tools (e.g., ANSYS Fluent) is recommended.

### Physical Principles
The calculator estimates the view factor, quantifying the radiation leaving one surface and striking another. Key components:
- **Emission Surface (Body A)**: Rays originate from Body A and aim at the target.
- **Target Surface (Body B)**: Receives rays from Body A.
- **Obstacles**: Optional, blocking rays and affecting the view factor.

### Features
#### Core Functionalities
- **Surface Loading**: Import Body A, Body B, and optional obstacles from STL files.
- **Calculation Modes**: Options include two-body interaction, body-environment interaction, and interactions with obstacles.
- **Ray Tracing**: Simulate with specified ray quantities.
- **Data Handling**: Save results for analysis and load saved data.

#### Visualization Options
- **3D Surface Visualization**: Display loaded surfaces.
- **Ray Path Visualization**: View ray paths and interactions.

### Installation and Execution
1. Save the MATLAB class as `ViewFactorCalculatorAppV5.m`.
2. Open MATLAB, navigate to the file's directory, and run:
   ```matlab
   app = ViewFactorCalculatorAppV5;
