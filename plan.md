Global Stability is the toughest part, I have structured this proposal to use a standard **Volterra-type Lyapunov function**. This means we won't be guessing the proof from scratch; we'll be adapting a known mathematical template.

---

# Project Proposal & Execution Blueprint

**To:** Instructor Asgher, Advanced Differential Equations 
**Topic:** *Survival of the Fittest: A Mathematical Study of Multi-Strain SIR Dynamics under Differential Vaccine Efficacy*
**Scope:** Local Stability, Global Stability, Numerical Simulation, and Sensitivity Analysis

---

## 1. Project Overview and Mathematical Model
This project investigates the competitive dynamics of two distinct viral strains (e.g., SARS-CoV-2 Delta vs. Omicron) competing for a single host population. We will extend the standard epidemiological framework to an **SVIR (Susceptible-Vaccinated-Infected-Recovered)** model to analyze how differential vaccine efficacy ($\epsilon$) dictates which strain achieves endemicity and which is driven to extinction (The Competitive Exclusion Principle).

**The Governing System of Non-Linear ODEs:**
$$\frac{dS}{dt} = \Lambda - (\beta_1 I_1 + \beta_2 I_2)S - (\nu + \mu)S$$
$$\frac{dV}{dt} = \nu S - (\sigma_1 \beta_1 I_1 + \sigma_2 \beta_2 I_2)V - \mu V$$
$$\frac{dI_1}{dt} = (\beta_1 S + \sigma_1 \beta_1 V)I_1 - (\gamma_1 + \mu)I_1$$
$$\frac{dI_2}{dt} = (\beta_2 S + \sigma_2 \beta_2 V)I_2 - (\gamma_2 + \mu)I_2$$
$$\frac{dR}{dt} = \gamma_1 I_1 + \gamma_2 I_2 - \mu R$$

*(Where $\sigma_n = 1 - \epsilon_n$, representing the "vaccine leakage" or immune escape of strain $n$.)*

---

## Phase 1: Local Stability Analysis
**High-Level Goal:** Determine the baseline mathematical conditions under which the disease naturally dies out or becomes an epidemic.

**Low-Level Execution (Step-by-Step):**
1.  **Calculate the Disease-Free Equilibrium (DFE):** 
    *   Set all derivatives to $0$. 
    *   Set $I_1 = 0$ and $I_2 = 0$.
    *   Solve algebraically for $S^*, V^*,$ and $R^*$ to yield the DFE state: $E_0$.
2.  **Construct the Jacobian Matrix ($J$):**
    *   Compute the partial derivatives of all 5 equations with respect to all 5 variables ($S, V, I_1, I_2, R$), resulting in a $5 \times 5$ matrix.
3.  **Evaluate $J$ at $E_0$:**
    *   Substitute the coordinates of the DFE into the Jacobian. This simplifies the matrix significantly, creating a block-triangular structure.
4.  **Derive the Basic Reproduction Number ($R_0$):**
    *   Extract the eigenvalues from the evaluated Jacobian.
    *   Alternatively, use the **Next-Generation Matrix Method** ($F \cdot V^{-1}$) on the infected compartments ($I_1, I_2$) to clearly define $R_{0,1}$ and $R_{0,2}$.
    *   *Proof objective:* Establish that if $\max(R_{0,1}, R_{0,2}) < 1$, the DFE is locally asymptotically stable.

---

## Phase 2: Global Stability Analysis
**High-Level Goal:** Prove mathematically that regardless of the initial number of infected people, the system will always eventually converge to the equilibrium state over infinite time ($t \to \infty$). 

**Low-Level Execution (Step-by-Step):**
1.  **Select the Lyapunov Function:**
    *   We will construct a standard Volterra-type Lyapunov function for the DFE. The general template to use is:
        $$L(S, V, I_1, I_2) = (S - S^* - S^* \ln \frac{S}{S^*}) + (V - V^* - V^* \ln \frac{V}{V^*}) + c_1 I_1 + c_2 I_2$$
2.  **Compute the Orbital Derivative ($\frac{dL}{dt}$):**
    *   Apply the chain rule: $\frac{dL}{dt} = \frac{\partial L}{\partial S}\frac{dS}{dt} + \frac{\partial L}{\partial V}\frac{dV}{dt} + \dots$
    *   Substitute the original ODEs into this derivative.
3.  **Algebraic Manipulation:**
    *   Group the terms and simplify using the equilibrium identities (e.g., $\Lambda = (\nu + \mu)S^*$).
    *   Apply algebraic inequalities (such as $1 - x + \ln x \le 0$ for $x > 0$) to prove that all terms in the derivative are strictly negative or zero.
4.  **Apply LaSalle's Invariance Principle:**
    *   Conclude the proof by stating that since $\frac{dL}{dt} \le 0$, the largest invariant set where $\frac{dL}{dt} = 0$ is the DFE itself, proving global asymptotic stability when $R_0 \le 1$.

---

## Phase 3: Computational Simulation
**High-Level Goal:** Numerically integrate the ODE system to visualize the population dynamics over a 365-day period, bridging the theoretical math with applied data.

**Low-Level Execution (Step-by-Step):**
1.  **Environment Setup:**
    *   Initialize a Python environment with `numpy`, `scipy.integrate`, and `matplotlib.pyplot`.
2.  **Parameterization:**
    *   Hardcode the parameter arrays. Use real-world estimates: $\beta_1 = 0.60$ (Delta), $\beta_2 = 1.16$ (Omicron), $\gamma_1 = \gamma_2 = 1/14$, $\nu = 0.01$.
    *   Set initial conditions: $S(0) = 0.948$, $V(0) = 0.05$, $I_1(0) = 0.001$, $I_2(0) = 0.001$, $R(0) = 0.0$.
3.  **Numerical Integration:**
    *   Define the ODE system as a Python function `def deriv(y, t, N, beta1, beta2, ...):`
    *   Use `scipy.integrate.odeint` or `solve_ivp` (RK45 method) to solve the array over `t = np.linspace(0, 365, 365)`.
4.  **Visualization:**
    *   Generate a multi-line time-series plot.
    *   Ensure academic formatting: grid lines, high DPI, distinct line styles (solid, dashed), and a clear legend. 

---

## Phase 4: Sensitivity Analysis
**High-Level Goal:** Demonstrate "originality and critical thinking" by systematically altering parameters to observe how vaccine evasion and transmission rates trigger competitive exclusion.

**Low-Level Execution (Step-by-Step):**
1.  **Scenario A (Baseline Competition):**
    *   Run the simulation with $\sigma_1 = 0.10$ (90% efficacy) and $\sigma_2 = 0.60$ (40% efficacy).
    *   *Expected output:* Strain 1 peaks early, but Strain 2 eventually dominates and drives Strain 1 to zero.
2.  **Scenario B (Vaccination Surge):**
    *   Increase the vaccination rate ($\nu$) from $0.01$ to $0.05$.
    *   *Expected output:* The peak of Strain 1 is flattened entirely, but Strain 2 still surges due to high immune escape ($\sigma_2$), proving that vaccines alone cannot stop a highly evasive variant.
3.  **Scenario C (Parameter Sweeping / Phase Portrait):**
    *   Write a loop to run the simulation across 20 different values of $\beta_2$ (ranging from $0.5$ to $1.5$).
    *   Plot the *final size* of the $I_2$ outbreak against the $\beta_2$ values to visually show the exact threshold where the variant takes over.

---

### Team Workload Distribution (Estimated 20 Cumulative Hours)
To execute this efficiently, divide the tasks as follows:

*   **Student 1 (The Analyst - 10 Hours):** 
    *   Executes Phase 1 (Local Stability).
    *   Executes Phase 2 (Global Stability). 
    *   Writes the mathematical proofs section of the final LaTeX report.
*   **Student 2 (The Programmer - 10 Hours):** 
    *   Executes Phase 3 (Simulation Code).
    *   Executes Phase 4 (Sensitivity Analysis scenarios).
    *   Formats the graphs, inserts them into the report, and writes the computational discussion.