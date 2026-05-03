import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Global constant for natural mortality rate (Assuming average lifespan of 70 years)
MU = 1 / (70 * 365)

def svir_model(y, t, Lambda, beta1, beta2, nu, mu, sigma1, sigma2, gamma1, gamma2):
    """
    SVIR (Susceptible-Vaccinated-Infected-Recovered) epidemiological model with two strains.
    
    Parameters:
    y : list
        Initial conditions [S, V, I1, I2, R]
    t : array
        Time points to solve the ODEs at
    Lambda : float
        Recruitment rate (days^-1)
    beta1, beta2 : float
        Transmission rates for Strain 1 and Strain 2 (days^-1)
    nu : float
        Vaccination rate (days^-1)
    mu : float
        Natural mortality rate (days^-1)
    sigma1, sigma2 : float
        Vaccine leakage/immune escape for Strain 1 and Strain 2 (dimensionless, 0 to 1)
    gamma1, gamma2 : float
        Recovery rates for Strain 1 and Strain 2 (days^-1)
        
    Returns:
    list
        Derivatives [dSdt, dVdt, dI1dt, dI2dt, dRdt]
    """
    S, V, I1, I2, R = y
    
    dSdt = Lambda - (beta1 * I1 + beta2 * I2) * S - (nu + mu) * S
    dVdt = nu * S - (sigma1 * beta1 * I1 + sigma2 * beta2 * I2) * V - mu * V
    dI1dt = (beta1 * S + sigma1 * beta1 * V) * I1 - (gamma1 + mu) * I1
    dI2dt = (beta2 * S + sigma2 * beta2 * V) * I2 - (gamma2 + mu) * I2
    dRdt = gamma1 * I1 + gamma2 * I2 - mu * R
    
    return [dSdt, dVdt, dI1dt, dI2dt, dRdt]

def run_simulation(beta1=0.60, beta2=1.16, nu=0.01, sigma1=0.10, sigma2=0.60, days=1095):
    # Fixed parameters
    mu = MU 
    Lambda = mu # Assuming constant population size N=1 => Lambda = mu * N
    gamma1 = 1 / 14
    gamma2 = 1 / 14
    
    # Initial conditions
    S0 = 0.948
    V0 = 0.05
    I10 = 0.001
    I20 = 0.001
    R0 = 0.0
    y0 = [S0, V0, I10, I20, R0]
    
    # Finer time resolution (10 points per day)
    t = np.linspace(0, days, days * 10)
    
    # Run ODE solver
    sol = odeint(svir_model, y0, t, args=(Lambda, beta1, beta2, nu, mu, sigma1, sigma2, gamma1, gamma2))
    
    return t, sol

def plot_simulation(t, sol, title, filename=None):
    S, V, I1, I2, R = sol.T
    
    plt.figure(figsize=(10, 6), dpi=300)
    plt.plot(t, S, label='Susceptible (S)', linestyle='-')
    plt.plot(t, V, label='Vaccinated (V)', linestyle='-')
    plt.plot(t, I1, label='Strain 1 Infected (Delta)', linestyle='--', color='red')
    plt.plot(t, I2, label='Strain 2 Infected (Omicron)', linestyle='-.', color='purple')
    plt.plot(t, R, label='Recovered (R)', linestyle='-')
    
    plt.xlabel('Days')
    plt.ylabel('Population Fraction')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    
    if filename:
        plt.savefig(filename, bbox_inches='tight')
        print(f"Saved plot to {filename}")
    else:
        plt.show()
    plt.close()

if __name__ == "__main__":
    # Scenario A: Baseline Competition
    # Explicitly setting days=1095 to align with the report's endemicity claim
    t, sol_A = run_simulation(sigma1=0.10, sigma2=0.60, nu=0.01, days=1095)
    plot_simulation(t, sol_A, "Scenario A: Baseline Competition", "scenario_A.png")
    
    # Scenario B: Vaccination Surge
    t, sol_B = run_simulation(sigma1=0.10, sigma2=0.60, nu=0.05, days=1095)
    plot_simulation(t, sol_B, "Scenario B: Vaccination Surge", "scenario_B.png")

    # Scenario C: Parameter Sweeping / Phase Portrait
    beta2_values = np.linspace(0.0, 1.5, 30)
    peak_I2_prevalence = []
    
    for b2 in beta2_values:
        t, sol = run_simulation(beta2=b2, days=500) # Run slightly longer to ensure equilibrium
        # Take the peak prevalence of I2
        max_I2 = np.max(sol[:, 3])
        peak_I2_prevalence.append(max_I2)
        
    plt.figure(figsize=(10, 6), dpi=300)
    plt.plot(beta2_values, peak_I2_prevalence, marker='o', linestyle='-', color='purple')
    plt.xlabel(r'Transmission Rate of Strain 2 ($\beta_2$)')
    plt.ylabel('Peak Prevalence of Strain 2')
    plt.title('Scenario C: Sensitivity of Strain 2 Peak Prevalence to Transmission Rate')
    plt.grid(True)
    plt.savefig("scenario_C.png", bbox_inches='tight')
    print("Saved plot to scenario_C.png")
