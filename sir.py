import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the SIR model
def sir_model(t, y, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]

# Parameters
beta = 0.8   # Infection rate
gamma = 0.1  # Recovery rate
N = 1000     # Total population

# Initial conditions
I0 = 1       # Initial infected
S0 = N - I0  # Initial susceptible
R0 = 0       # Initial recovered

Rr0 = beta / gamma 
print(Rr0)

y0 = [S0 / N, I0 / N, R0 / N]  # Normalize to fractions

# Time span
t_span = (0, 160)  # Simulate for 160 days
t_eval = np.linspace(*t_span, 1000)  # Time points for evaluation

# Solve the system
solution = solve_ivp(sir_model, t_span, y0, args=(beta, gamma), t_eval=t_eval)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[0], label="Susceptible", color='blue')
plt.plot(solution.t, solution.y[1], label="Infected", color='red')
plt.plot(solution.t, solution.y[2], label="Recovered", color='green')
plt.xlabel("Time (days)")
plt.ylabel("Fraction of Population")
plt.legend()
plt.title("SIR Model - Kermack-McKendrick")
plt.grid()
plt.show()
