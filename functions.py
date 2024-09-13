def schwarzschild_orbit(a, ecc, N=10**7): 
    # a -> semi-major axis (in black hole mass units)
    # ecc -> eccentricity
    # N -> number of steps

    from scipy import integrate
    import numpy as np

    def V2(r,L):
        return (1-2/r)*(1+L**2/r**2)

    L=np.sqrt(a*(1-ecc**2))
    E2= V2(a*(1-ecc),L)

    def evol(t,x):
        theta, r, dotr = x
        fphi = L/r**2
        fr = dotr
        fdotr = -1/(r**2)+(L**2)*(1-3/r)/r**3
        return fphi, fr, fdotr
    
    r0=a              # initial distance
    vr20 = E2-V2(r0,L)  # initial squared radial velocity
    if vr20 <= 0:
        raise ValueError("Particle too close to black hole, try increasing 'a' or decreasing 'ecc'")
    vr0=np.sqrt(vr20)      # initial radial velocity
    phi0=0          # initial angular coordinate

    T = (200*np.pi)*r0**2/L  # simulation time (around 5 orbits)
    
    def sol(T, phi0, r0, vr0, N):
        orbit = integrate.solve_ivp(evol,(0,T),(phi0,r0, vr0), t_eval=np.linspace(0,T,N))
        phi, r, vr = np.array(orbit.y)
        return phi, r
    
    return sol(T, phi0, r0, vr0, N)

def schwarzschild_plot_orbit(a, ecc, N=10**7):
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Gerar dados da órbita
    phi, r = schwarzschild_orbit(a, ecc, N)

    phi_p = phi[phi < 10]
    r_p = r[phi < 10]
        
    
    # Configurações do gráfico
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ax.plot(phi_p, r_p, color="r", linewidth=1.5, linestyle='-')
    # Adicionar o círculo representando o buraco negro
    black_hole = plt.Circle((0, 0), 1, transform=ax.transData._b, color='black', alpha=1)
    ax.add_artist(black_hole)
    
    # Ajustar limites de r
    ax.set_ylim(0, max(r)*1.1)
    
    # Remover legendas radiais e angulares
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    # Título
    ax.set_title('Trajetória da Partícula em torno do Buraco Negro', va='bottom', fontsize=14)
    
    # Customização da grid e fundo
    ax.grid(True, color="gray", linestyle="--", linewidth=0.5)
    ax.set_facecolor("#1a1a1a")  # Fundo mais claro
    
    plt.show()

    return 'Plot displayed successfully'

def perihelion(a, ecc, N=10**7):

    from scipy.signal import argrelextrema
    import numpy as np
    from sklearn import linear_model

    phi, r = schwarzschild_orbit(a, ecc)

    r_perid=argrelextrema(r,np.less)[0]

    phi_per = [phi[r_perid[i]] - 2*np.pi*i for i in range(len(r_perid))]

    
    reg = linear_model.LinearRegression()
    phi_perfit = np.arange(len(phi_per)).reshape(-1, 1)
    reg.fit(phi_perfit,phi_per)

    z=reg.coef_[0]

    return z


