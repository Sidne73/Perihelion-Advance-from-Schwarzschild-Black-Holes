# Black Hole Orbital Simulation

This project simulates the motion of a particle around a static black hole (Schwarzschild black hole) and plots its trajectory. The simulation is done using a 4th-order Runge-Kutta method to calculate the orbit from perihelion (the point of closest approach) to aphelion (the farthest point). Additionally, it calculates the perihelion precession due to general relativistic effects.

## Overview

### Functions

- **`V2(r, L)`**:  
  Calculates the effective potential for a given radial distance `r` and angular momentum `L`.

- **`elliptic(L)`**:  
  Returns the minimum energy squared (`E²_min`) for a closed orbit with a given angular momentum `L`.

- **`r_lim(L, p)`**:  
  Finds the radial limits (`u_min`, `u_max`), where `u = 1/r`, that bound the orbit for a given angular momentum `L` and parameter `p`. The parameter `p` is a measure of the orbits energy, when `p=0` we have `E=E_min` and the orbit will be circular, and when `p=1` we have `E=1`, the orbit will have a return point at infinity, resembling a parabola. Thus `p` is closely related to the orbits' eccentricity.

- **`schwarzschild_orbit(L, p, N)`**:  
  Simulates the orbit of a particle in Schwarzschild spacetime using the 4th-order Runge-Kutta method, returning the angular position `phi` and radial distance `r`.

- **`perihelion(L, p)`**:  
  Calculates the perihelion precession of the orbit, returning the angular shift due to relativistic effects.

- **`schwarzschild_plot_orbit(L, p, N)`**:  
  Plots the orbit of the particle around the black hole using polar coordinates. The orbit is extended and reflected to simulate multiple revolutions.

## Requirements

Make sure to install the required libraries:

```bash
pip install numpy scipy sympy matplotlib
```

## Usage

The main script can be run inside a Jupyter notebook to demonstrate the functions and visualize the results. Below is an example of how to use the functions in your `example.ipynb`.

### 1. Calculating the Perihelion Precession

You can calculate the perihelion precession for a particle orbiting the black hole:

```python
L = 50.0  # Angular momentum
p = 0.5  # Parameter related to eccentricity

precession = perihelion(L, p)
print(f"Perihelion Precession: {precession} radians/orbit")
```


### 2. Plotting the Orbit
You can plot the orbit of a particle around the black hole using the following command:

```python
L = 4.0  # Angular momentum
p = 0.1  # Parameter related to eccentricity

schwarzschild_plot_orbit(L, p)
```

This will display a polar plot of the orbit, showing the trajectory of the particle around the black hole. The black hole is represented as a dark circle in the center, and the orbit is extended over multiple revolutions for visualization.

### 3. Modifying the Simulation
The `schwarzschild_orbit(L, p, N)` function allows you to simulate orbits with a different resolution by changing the value of `N`. Increasing `N` will increase the accuracy of the simulation but may take more time to compute:

```python
L = 4.0  # Angular momentum
p = 0.1  # Parameter related to eccentricity
N = 10**6  # Higher resolution

phi, r = schwarzschild_orbit(L, p, N)
```
## Project Structure
```bash
.
├── example.ipynb    # Jupyter notebook with examples of usage
├── functions.py    # Python file with the orbital simulation functions
└── README.md        # This documentation file
```
## Running the Project

1. Open the `example.ipynb` notebook.

2. Execute the cells to run the examples provided for calculating perihelion precession and plotting orbits.

3. Modify the parameters `L`, `p`, and `N` as needed to explore different orbital configurations. Due to the way this program works, simulations for large values of `L` (`L>1000`) or low values of `p` (this may depend on `L`) may lead to inaccurate results.


## 🤝 Contributions and Contact

If you'd like to contribute to this project, feel free to contact me at:  
📧 **sidney1395271@gmail.com**

## 🔮 Future Work

This simulation is very limited in its capabilities, there must be a way to simulate the trajectories with high precision with a wider range of parameters than this provides.

