import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

def plot_graph():
    # Define x values
    L = 85 * 10**-10
    x_values_1 = np.linspace(-L/2, L/2, 1000)
    x_values_2 = np.linspace(-L, -L/2, 500)
    x_values_3 = np.linspace(L/2, L, 500)
    
    # Calculate ground state energy for the structure
    E_finite = 0.01743222976650539 * 1.6022 *10**-19

    fig, ax = plt.subplots()
    
    # Define y values
    m = 0.038 * 9.109 * 10**-31
    V_0 = 0.15403769200000003 * 1.6022 *10**-19
    h_bar_2 = (6.62607015 * 10**-34 / (2 * math.pi))**2 
    E_infinite = 2.194381046748574 *10**-20
    a = math.sqrt(abs((2 * m * (V_0 - E_infinite) / h_bar_2)))
    k = math.sqrt((2*m*E_infinite)/h_bar_2)
    y_cosine = (np.cos(k * (x_values_1))) * (E_finite) + E_finite
    y_exponential_positive = (E_finite * np.exp(a * x_values_2)) * np.exp(a * L/2)
    y_exponential_negative = (E_finite * np.exp(-a * x_values_3)) * np.exp(a * L/2)
    #ax.hlines(y= E_finite, xmin=-L/2, xmax=L/2, color='green', label='Ground state energy')
    ax.hlines(y= V_0, xmin=-L, xmax=-L/2, color='red', label='Valence band')
    ax.hlines(y= 0, xmin=-L/2, xmax=L/2, color='red')
    ax.hlines(y= V_0, xmin=L/2, xmax=L, color='red')
    ax.vlines(x=-L/2, ymin=0, ymax=V_0, color='red')
    ax.vlines(x=L/2, ymin=0, ymax=V_0, color='red')
    ax.text(0, 0.75*10**(-20), f'Ground state energy: {E_finite:5g} J', ha='center', va='center', fontsize=8, color='black')
    
    # Plot the graphs
    plt.plot(x_values_1, y_cosine, color = 'green', label = 'Wavefunction')
    plt.plot(x_values_2, y_exponential_positive, color ='green')
    plt.plot(x_values_3, y_exponential_negative, color = 'green')
    plt.xlabel('L (m^-1)')
    plt.ylabel('Energy (J)')
    plt.title('Finite potential well even solution for electrons')
    plt.legend()
    plt.grid(True)
    plt.show()

plot_graph()

