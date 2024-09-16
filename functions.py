# Function that calculates the effective potential for a given radial distance (r) and angular momentum (L)
def V2(r, L):
    return (1 - 2/r) * (1 + L**2 / r**2)

# Function to determine the minimum energy squared (E²_min) for an orbit with a given angular momentum (L)
def elliptic(L):
    import numpy as np
    # Calculate the radius of the inner circular orbit (r_circ1)
    r_circ1 = L**2 * (1 + np.sqrt(1 - 12 / L**2)) / 2
    # Calculate the minimum energy squared (E²_min) at this orbit
    E2_min = V2(r_circ1, L)
    return E2_min

# Function that determines the minimum and maximum radial limits (u_min, u_max) of the orbit for a given L and p
def r_lim(L, p):
    import sympy
    from sympy import symbols, solve, N
    from sympy.functions.elementary.complexes import re

    # Compute the total energy squared (E²) for the given parameters
    E2 = (1 - p) * elliptic(L) + p

    # Define symbolic variables for energy (E²), angular momentum (L), and the inverse of the radial distance (u = 1/r)
    E2_sym = symbols('E2_sym')
    L_sym = symbols('L_sym')
    u_sym = symbols('u_sym')

    # Solve the equation for the radial limits (u) based on the effective potential
    u_lim = solve(E2_sym - (1 - 2 * u_sym) * (1 + L_sym**2 * u_sym**2), u_sym)

    # Extract the real parts of the solutions to obtain the minimum (u_min) and maximum (u_max) values of u
    u_max = re(N(u_lim[1].subs(E2_sym, E2).subs(L_sym, L)))
    u_min = re(N(u_lim[0].subs(E2_sym, E2).subs(L_sym, L)))

    # Return the radial limits as floating point numbers
    return float(u_min), float(u_max)

# Main function to simulate the orbit of a particle in the Schwarzschild spacetime for given angular momentum (L) and parameter (p)
def schwarzschild_orbit(L, p, N=10**5):
    from scipy import integrate
    import numpy as np
    import matplotlib.pyplot as plt

    # Compute the total energy squared (E²) for the given parameters
    E2 = (1 - p) * elliptic(L) + p

    # Get the minimum and maximum radial limits (u_min, u_max)
    u_min, u_max = r_lim(L, p)

    # Function to compute the first derivative of u with respect to φ (angular coordinate)
    def dudphi(u):
        return np.sqrt(E2 / L**2 - (1 - 2 * u) * (u**2 + 1 / L**2))
    
    # Helper function for root-finding to check when the orbit reaches u_min and u_max
    def dudphi2(u):
        return E2 / L**2 - (1 - 2 * u) * (u**2 + 1 / L**2)
    
    # Ensure that u_min corresponds to a valid point where the derivative is near zero
    while dudphi2(u_min) > 0:
        u_min -= u_min / 10000000000
    
    # Ensure that u_max corresponds to a valid point where the derivative is near zero
    while dudphi2(u_max) > 0:
        u_max += u_max / 10000000000
    
    # Ensuring the derivative is not precisely zero because this would result in a static solution
    while dudphi2(u_min) < u_min * 10**(-14):
        u_min += u_min * 10**(-14)
    
    while dudphi2(u_max) < u_max * 10**(-14):
        u_max -= u_max * 10**(14)

    # Function to evolve the orbit using the 4th order Runge-Kutta method
    def evol(u_0, u_1, N):
        # Step size for φ (angular coordinate) is π/N
        step = (np.pi) / N
        u = np.zeros(10 * N)  # Array to store values of u (1/r)
        u[0] = u_0  # Initial value of u
        i = 0
        # Evolve u from u_0 (perihelion) to u_1 (aphelion)
        while u[i] < u_1:
            # Compute the four stages of Runge-Kutta method
            k1 = dudphi(u[i])
            k2 = dudphi(u[i] + k1 * step / 2)
            k3 = dudphi(u[i] + k2 * step / 2)
            k4 = dudphi(u[i] + k3 * step)
            # Calculate the slope as a weighted average of the four stages
            slope = (k1 + 2 * k2 + 2 * k3 + k4) / 6
            u[i+1] = u[i] + step * slope  # Update the value of u
            i += 1
        # Compute the corresponding φ values and return φ and u
        phi = np.arange(0, i * step, step)
        return phi, u[:i]
    
    # Evolve the orbit from perihelion (u_min) to aphelion (u_max)
    phi, u = evol(u_min, u_max, N)
    
    # Convert u (1/r) to r (radial distance) for the final result
    r = [1 / x for x in u]
    return phi, r
        
# Function to calculate the perihelion precession for a given angular momentum (L) and parameter (p)
def perihelion(L, p):
    import numpy as np
    # Compute the orbit data (angular position φ and radial distance r)
    phi, r = schwarzschild_orbit(L, p)
    # Return the perihelion precession: 2 * (final angle - π)
    return 2 * (phi[-1] - np.pi)

# Function to plot the Schwarzschild orbit using a polar plot
def schwarzschild_plot_orbit(L, p, N=10**5):
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Generate the orbit data (angular position φ and radial distance r)
    phi, r = schwarzschild_orbit(L, p, N)
    step = (np.pi) / N  # Step size for φ

    # Create a range of additional angular positions to extend the plot
    phirange = np.arange(1, len(phi))
    phi = np.append(phi, phi[-1] + phirange * step)  # Extend φ array
    r = np.concatenate((r, r[-2::-1]))  # Reflect the r array to simulate one full orbit
    
    # Number of additional orbits to plot (for visualization)
    n_orbits = 5
    phirange = np.arange(n_orbits * len(phi))  # Create angular positions for multiple orbits
    phi = np.append(phi, phi[-1] + phirange * step)  # Extend φ to cover multiple orbits
    r2 = r
    for i in range(n_orbits):
        # Concatenate the reflected r values to simulate multiple orbits
        r = np.concatenate((r, r2))
    
    # Plot configuration
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})  # Create a polar plot
    ax.plot(phi, r, color="r", linewidth=1.5, linestyle='-')  # Plot the orbit with red lines
    
    # Add a circle representing the black hole at the center (radius 1)
    black_hole = plt.Circle((0, 0), 1, transform=ax.transData._b, color='black', alpha=1)
    ax.add_artist(black_hole)  # Add the black hole to the plot
    
    # Adjust the radial limit of the plot
    ax.set_ylim(0, max(r) * 1.1)  # Set radial limits to be slightly larger than max(r)
    
    # Remove the radial and angular labels for a cleaner plot
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    # Set the title of the plot
    ax.set_title('Orbit Plot', va='bottom', fontsize=14)
    
    # Customize the grid and background
    ax.grid(True, color="gray", linestyle="--", linewidth=0.5)  # Add a gray dashed grid
    ax.set_facecolor("#1a1a1a")  # Set the background color to a slightly lighter gray
    
    # Display the plot
    plt.show()

    return 'Plot displayed successfully'


