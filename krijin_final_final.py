import math
import pandas as pd
import matplotlib.pyplot as plt

# Read data from files
krijin_table1 = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Krijin Table 1.xlsx', index_col='Material')
krijin_table2 = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Krijin Table 2.xlsx', index_col='Material')
ingaas_table = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\InGaAs Table 1.xlsx', index_col='Structure')
data_001 = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Data_001.xlsx', index_col='Material')
data_111 = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Data_111.xlsx', index_col='Material')

def calculate_bandgap(substrate, epilayer, orientation):
    """
    Calculate the bandgap of the epilayer material on a given substrate for a specific crystal orientation.

    Parameters:
        substrate (tuple): Material parameters of the substrate.
        epilayer (tuple): Material parameters of the epilayer.
        orientation (str): Crystal orientation ('001' or '111').

    Returns:
        tuple: Valence band edge of the epilayer, conduction band edge of the epilayer,
               valence band edge of the substrate, conduction band edge of the substrate,
               energy gap between the valence bands of the epilayer and substrate.
    """
    # Material parameters extraction
    (a0_substrate, _, _, _, Evav_substrate,
     Δ0_substrate, Eg_Γ_substrate, _, _, _,
     _, _, _, _, _) = substrate

    (a_epilayer, c11_epilayer, c12_epilayer, c44_epilayer, Evav_epilayer,
     Δ0_epilayer, Eg_Γ_epilayer, _, _, av_epilayer,
     ac_Γ_epilayer, b_epilayer, d_epilayer, _, _) = epilayer

    # Equations 3a and 3b
    D_001 = 2 * c12_epilayer / c11_epilayer
    D_111 = 2 * (c11_epilayer + 2 * c12_epilayer - 2 * c44_epilayer) / (c11_epilayer + 2 * c12_epilayer + 4 * c44_epilayer)

    # Determine lattice constants based on orientation using Equations 2a and 2b
    if orientation == '001':
        a_parallel_substrate = a0_substrate
        a_perpendicular_substrate = a_epilayer * (1 - D_001 * (a_parallel_substrate / a_epilayer - 1))
    elif orientation == '111':
        a_parallel_substrate = a0_substrate
        a_perpendicular_substrate = a_epilayer * (1 - D_111 * (a_parallel_substrate / a_epilayer - 1))

    # Determine strain on interface (Equation 1)
    ε_parallel = a_parallel_substrate / a_epilayer - 1
    ε_perpendicular = a_perpendicular_substrate / a_epilayer - 1

    # Hydrostatic strain contributions for valence and conduction bands (Equations 4a and 4b)
    Δ_Evav_hydrostatic = av_epilayer * (2 * ε_parallel + ε_perpendicular)
    Δ_Ec_hydrostatic = ac_Γ_epilayer * (2 * ε_parallel + ε_perpendicular)

    # Shear strain contributions (Equations 6a and 6b)
    if orientation == '001':
        δ_E_sh = 2 * b_epilayer * (ε_perpendicular - ε_parallel)
    elif orientation == '111':
        δ_E_sh = (2 / 3) * math.sqrt(3) * d_epilayer * (ε_perpendicular - ε_parallel)

    # Equations 5a, and 5b (5c not used)
    Δ_Ehh_sh = -0.5 * δ_E_sh
    Δ_Elh_sh = -0.5 * Δ0_epilayer + 0.25 * δ_E_sh + \
                0.5 * math.sqrt(Δ0_epilayer**2 + Δ0_epilayer * δ_E_sh + 9 / 4 * δ_E_sh**2)

    # Calculate valence and conduction band edges for the epilayer using Equations 7a and 7b
    Ev_epilayer = Evav_epilayer + Δ0_epilayer / 3 + Δ_Evav_hydrostatic + max(Δ_Ehh_sh, Δ_Elh_sh)
    Ec_epilayer = Evav_epilayer + Δ0_epilayer / 3 + Eg_Γ_epilayer + Δ_Ec_hydrostatic

    # Calculate valence and conduction band edges for the substrate
    Ev_substrate = Evav_substrate + Δ0_substrate / 3
    Ec_substrate = Evav_substrate + Δ0_substrate / 3 + Eg_Γ_substrate

    # Calculate the energy gap between the valence band of the epilayer and the substrate
    ΔE_valence_gap = Ev_epilayer - Ev_substrate

    return Ev_epilayer, Ec_epilayer, Ev_substrate, Ec_substrate, ΔE_valence_gap, D_001, a_parallel_substrate, a_perpendicular_substrate, ε_parallel, ε_perpendicular, Δ_Evav_hydrostatic, Δ_Ec_hydrostatic, δ_E_sh, Δ_Ehh_sh, Δ_Elh_sh

def plot_binary(ax, ΔE_valence_gap, substrate_material, epilayer_material, Ev_substrate, Ev_epilayer, Ec_substrate, Ec_epilayer, orientation):
    """
    Plot the energy levels for binary epilayer-substrate combinations.

    Parameters:
        ax (matplotlib.axes.Axes): The subplot axes.
        ΔE_valence_gap (float): Energy gap between the valence bands of the epilayer and substrate.
        substrate_material (str): Name of the substrate material.
        epilayer_material (str): Name of the epilayer material.
        Ev_substrate (float): Valence band edge of the substrate.
        Ev_epilayer (float): Valence band edge of the epilayer.
        Ec_substrate (float): Conduction band edge of the substrate.
        Ec_epilayer (float): Conduction band edge of the epilayer.
        orientation (str): Crystal orientation ('001' or '111').
    """
    ax.bar([0], [Ec_substrate - Ev_substrate], bottom=Ev_substrate, width=0.4, color='blue', label=substrate_material)
    ax.bar([0.4], [Ec_epilayer - Ev_epilayer], bottom=Ev_epilayer, width=0.4, color='red', label=epilayer_material)
    ax.set_ylabel('Energy (eV)')
    ax.set_xlabel('Materials')
    ax.set_title(f'Energy Levels of {substrate_material} Substrate and {epilayer_material} Epilayer ({orientation} Orientation)')
    ax.set_xticks([])  # Remove x-axis ticks
    ax.legend()

    # Annotate the energy gap between the valence band of the epilayer and the valence band of the substrate
    ΔE_annotation = ΔE_valence_gap
    arrow_start = (0.4, Ev_substrate) 
    arrow_end = (0.4, Ev_epilayer)

    ax.text(0.51, (Ev_substrate + Ev_epilayer) / 2, 
            f'Energy Gap: {ΔE_annotation:.2f} eV', 
            ha='center', va='center', fontsize=8, color='black')
    ax.annotate(f'',
                xy=arrow_start, xycoords='data',
                xytext=arrow_end, textcoords='data',
                arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"),
                fontsize=8, ha='center', va='top')

    # Annotate the valence band and the conduction band
    ax.text(-0.2, Ec_substrate, 'Ec', ha='right', va='center',
            fontsize=8, color='black')
    ax.text(-0.2, Ev_substrate, 'Ev', ha='right', va='bottom',
            fontsize=8, color='black')
    ax.text(0.6, Ec_epilayer, 'Ec', ha='left', va='center',
            fontsize=8, color='black')
    ax.text(0.6, Ev_epilayer, 'Ev', ha='left', va='center',
            fontsize=8, color='black')

def process_data_binary(orientations, epilayer_material, substrate_material):
    """
    Process data for binary epilayer-substrate combinations.

    Parameters:
        orientation (str): Crystal orientation ('001' or '111').

    Returns:
        tuple: The matplotlib figure and axes.
    """
    # Get material parameters for the specified epilayer and substrate
    epilayer_data = krijin_table1.loc[epilayer_material].values
    substrate_data = krijin_table1.loc[substrate_material].values

    # Calculate bandgap energies for orientations
    Ev_epilayer, Ec_epilayer, Ev_substrate, Ec_substrate, ΔE_valence_gap, D_001, a_parallel_substrate, a_perpendicular_substrate, ε_parallel, ε_perpendicular, Δ_Evav_hydrostatic, Δ_Ec_hydrostatic, δ_E_sh, Δ_Ehh_sh, Δ_Elh_sh = calculate_bandgap(substrate_data, epilayer_data, orientation=orientations)

    # Print the results
    print(f"Valence Band Edge of Substrate (Ev_substrate): {Ev_substrate} eV")
    print(f"Conduction Band Edge of Substrate (Ec_substrate): {Ec_substrate} eV")

    print(f"\nFor ({orientations}) orientation:")
    print(f"Valence Band Edge (Ev): {Ev_epilayer} eV")
    print(f"Conduction Band Edge (Ec): {Ec_epilayer} eV")
    print(f"Energy Gap between Valence Band of Epilayer and Substrate (ΔE_valence_gap): {ΔE_valence_gap} eV")

    if orientations == '001':
        data_001.loc[epilayer_material, 'Ev'] = Ev_epilayer
        data_001.loc[epilayer_material, 'Ec'] = Ec_epilayer
        data_001.loc[epilayer_material, 'D'] = D_001

        data_001.loc[epilayer_material, 'a_perpendicular'] = a_perpendicular_substrate
        data_001.loc[epilayer_material, 'ε_parallel'] = ε_parallel
        data_001.loc[epilayer_material, 'ε_perpendicular'] = ε_perpendicular
        data_001.loc[epilayer_material, 'Δ_Evav_hydrostatic'] = Δ_Evav_hydrostatic
        data_001.loc[epilayer_material, 'Δ_Ec_hydrostatic'] = Δ_Ec_hydrostatic
        data_001.loc[epilayer_material, 'δ_E_sh'] = δ_E_sh
        data_001.loc[epilayer_material, 'Δ_Ehh_sh'] = Δ_Ehh_sh
        data_001.loc[epilayer_material, 'Δ_Elh_sh'] = Δ_Elh_sh

        data_001.loc[substrate_material, 'Ev'] = Ev_substrate
        data_001.loc[substrate_material, 'Ec'] = Ec_substrate
        data_001.loc[substrate_material, 'a_parallel'] = a_parallel_substrate
        data_001.to_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Data_001.xlsx')
    
    elif orientations == '111':
        data_111.loc[epilayer_material, 'Ev'] = Ev_epilayer
        data_111.loc[epilayer_material, 'Ec'] = Ec_epilayer
        data_111.loc[epilayer_material, 'D_001'] = D_001
        data_111.loc[epilayer_material, 'a_perpendicular_substrate'] = a_perpendicular_substrate
        data_111.loc[epilayer_material, 'ε_parallel'] = ε_parallel
        data_111.loc[epilayer_material, 'ε_perpendicular'] = ε_perpendicular
        data_111.loc[epilayer_material, 'Δ_Evav_hydrostatic'] = Δ_Evav_hydrostatic
        data_111.loc[epilayer_material, 'Δ_Ec_hydrostatic'] = Δ_Ec_hydrostatic
        data_111.loc[epilayer_material, 'δ_E_sh'] = δ_E_sh
        data_111.loc[epilayer_material, 'Δ_Ehh_sh'] = Δ_Ehh_sh
        data_111.loc[epilayer_material, 'Δ_Elh_sh'] = Δ_Elh_sh

        data_111.loc[substrate_material, 'Ev'] = Ev_substrate
        data_111.loc[substrate_material, 'Ec'] = Ec_substrate
        data_111.loc[substrate_material, 'a_parallel'] = a_parallel_substrate
        data_111.to_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\data_111.xlsx')

    
    # Plot the results for orientations
    fig, ax = plt.subplots()
    plot_binary(ax, ΔE_valence_gap, substrate_material, epilayer_material, Ev_substrate, Ev_epilayer, Ec_substrate, Ec_epilayer, orientation=orientations)

    plt.show()

    return fig, ax

def calculate_T_values_ternary(structure, GaAs, InAs, GaInAs, material_parameter):
    """
    Calculate T values for ternary alloys.

    Parameters:
        structure (str): Structure name.
        GaAs (float): Material parameter for GaAs.
        InAs (float): Material parameter for InAs.
        GaInAs (float): Material parameter for GaInAs.
        material_parameter (str): Material parameter to calculate (Eg(Γ) or Δ0).

    Returns:
        numpy.ndarray: T values for the specified material parameter.
    """
    x_values = ingaas_table.loc['Ga mole fraction (x)', structure]
    if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0':
        T_ABC_values_ternary = x_values * GaAs + (1 - x_values) * InAs + x_values * (x_values - 1) * GaInAs
    elif material_parameter == 'Ev' or material_parameter == 'Ec':
        T_ABC_values_ternary = x_values * GaAs + (1 - x_values) * InAs

    return T_ABC_values_ternary

def process_data_ternary(material_parameter):
    """
    Process data for ternary alloys and return T values.

    Parameters:
        material_parameter (str): Material parameter to calculate (Ev, Ec, Eg(Γ)).

    Returns:
        list: T values for each structure.
    """
    # Get the values from Krijin Table 1 and 2
    GaAs = krijin_table1.loc["GaAs", material_parameter]
    InAs = krijin_table1.loc["InAs", material_parameter]
    if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0':
        GaInAs = krijin_table2.loc["GaxIn1-xAs", f'C({material_parameter})']

    # Initialize lists to store T_ABC values for each structure
    T_ABC_values_list = []

    # Calculate T_ABC for each structure and display the results
    structures = ingaas_table.columns[0:].tolist()
    for structure in structures:
        if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0':
            T_ABC_values = calculate_T_values_ternary(structure, GaAs, InAs, GaInAs, material_parameter)
        elif material_parameter == 'Ev' or material_parameter == 'Ec':
            GaInAs = None
            T_ABC_values = calculate_T_values_ternary(structure, GaAs, InAs, GaInAs, material_parameter)

        print(f'T_ABC_{structure}: {material_parameter}: {T_ABC_values} eV')
        T_ABC_values_list.append(T_ABC_values)

    return T_ABC_values_list

def calculate_T_values_quaternary(GaAs, InAs, GaInAs, GaP, GaInP, GaPAs, InP, InPAs, material_parameter):
    """
    Calculate T values for quaternary alloys.

    Parameters:
        GaAs (float): Material parameter for GaAs.
        InAs (float): Material parameter for InAs.
        GaInAs (float): Material parameter for GaInAs.
        GaP (float): Material parameter for GaP.
        GaInP (float): Material parameter for GaInP.
        GaPAs (float): Material parameter for GaPAs.
        InP (float): Material parameter for InP.
        InPAs (float): Material parameter for InPAs.
        material_parameter (str): Material parameter to calculate (Eg(Γ) or Δ0).

    Returns:
        tuple: T values for different compositions.
    """
    if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0':
        T_ABC_values_quaternary = 0.5 * GaAs + (1 - 0.5) * InAs + 0.5 * (0.5 - 1) * GaInAs
        T_ABD_values_quaternary = 0.5 * GaP + (1 - 0.5) * InP + 0.5 * (0.5 - 1) * GaInP
        T_ACD_values_quaternary = 0.5 * GaAs + (1 - 0.5) * GaP + 0.5 * (0.5 - 1) * GaPAs
        T_BCD_values_quaternary = 0.5 * InAs + (1 - 0.5) * InP + 0.5 * (0.5 - 1) * InPAs
    elif material_parameter == 'Ev' or material_parameter == 'Ec':
        T_ABC_values_quaternary = 0.5 * GaAs + (1 - 0.5) * InAs 
        T_ABD_values_quaternary = 0.5 * GaP + (1 - 0.5) * InP 
        T_ACD_values_quaternary = 0.5 * GaAs + (1 - 0.5) * GaP 
        T_BCD_values_quaternary = 0.5 * InAs + (1 - 0.5) * InP 

    return T_ABC_values_quaternary, T_ABD_values_quaternary, T_ACD_values_quaternary, T_BCD_values_quaternary

def calculate_Q_values(T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values):
    """
    Calculate Q values.

    Parameters:
        T_ABC_values (numpy.ndarray): T values for ABC composition.
        T_ABD_values (numpy.ndarray): T values for ABD composition.
        T_ACD_values (numpy.ndarray): T values for ACD composition.
        T_BCD_values (numpy.ndarray): T values for BCD composition.

    Returns:
        numpy.ndarray: Q values.
    """
    # Calculate Q(x, y) using Equation (10)
    Q_values = (
        (0.14 * (1 - 0.14) * (0.3 * T_ABC_values + (1 - 0.3) * T_ABD_values)) / (0.14 * (1 - 0.14) + 0.3 * (1 - 0.3)) +
        (0.3 * (1 - 0.3) * (0.14 * T_ACD_values + (1 - 0.14) * T_BCD_values)) / (0.14 * (1 - 0.14) + 0.3 * (1 - 0.3))
    )

    return Q_values

def process_data_quaternary(material_parameter):
    """
    Process data for quaternary alloys and return Q values.

    Parameters:
        material_parameter (str): Material parameter to calculate (Ev, Ec, Eg(Γ)).

    Returns:
        numpy.ndarray: Q values.
    """
    # Get material parameters
    GaAs = krijin_table1.loc["GaAs", material_parameter]
    InAs = krijin_table1.loc["InAs", material_parameter]
    GaP = krijin_table1.loc["GaP", material_parameter]
    InP = krijin_table1.loc["InP", material_parameter]
    if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0':
        GaInAs = krijin_table2.loc["GaxIn1-xAs", f'C({material_parameter})']
        GaInP = krijin_table2.loc["GaxIn1-xP", f'C({material_parameter})']        
        GaPAs  = krijin_table2.loc["GaPxAs1-x", f'C({material_parameter})'] 
        InPAs = krijin_table2.loc["InPxAs1-x", f'C({material_parameter})'] 
        T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values = calculate_T_values_quaternary(GaAs, InAs, GaInAs, GaP, GaInP, GaPAs, InP, InPAs, material_parameter)
    elif material_parameter == 'Ev' or material_parameter == 'Ec':
        GaInAs = GaInP = GaPAs = InPAs = None
        T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values = calculate_T_values_quaternary(GaAs, InAs, GaInAs, GaP, GaInP, GaPAs, InP, InPAs, material_parameter)
    # Calculate Q values for quaternary alloys
    Q_xy_values = calculate_Q_values(T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values)
    # Display the results
    print(f'Q(x,y): {material_parameter}: {Q_xy_values} eV ')

    return Q_xy_values
    
def plot_GaInAsP_InP(ax):
    """
    Plot the energy levels for GaInAsP/InP structures.

    Parameters:
        ax (matplotlib.axes.Axes): The subplot axes.
    """
    # Get the energy levels of InP
    InP_Ec = krijin_table1.loc["InP", 'Ec']
    InP_Ev = krijin_table1.loc["InP", 'Ev']

    # Get the energy levels of GaInAsP
    GaInAsP_Ec = process_data_quaternary('Ec')
    GaInAsP_Ev= process_data_quaternary('Ev')

    # Get T_ABC values for all structures for each material parameter
    structures = ingaas_table.columns[0:].tolist()
    GaInAs_Ec_list = process_data_ternary('Ec')
    GaInAs_Ev_list = process_data_ternary('Ev')

    # Iterate over structures
    for i, structure in enumerate(structures):
        GaInAs_Ec = GaInAs_Ec_list[i]
        GaInAs_Ev = GaInAs_Ev_list[i]

        # Create a new subplot for each structure
        fig, ax = plt.subplots()

        ax.bar([0], [InP_Ec - InP_Ev], bottom=InP_Ev, width=0.4, color='blue', label='InP')
        ax.bar([0.4], [GaInAsP_Ec - GaInAsP_Ev], bottom=GaInAsP_Ev, width=0.4, color='red', label='GaInAsP')
        ax.bar([0.8], [GaInAs_Ec - GaInAs_Ev], bottom=GaInAs_Ev, width=0.4, color='green', label='GaInAs')
        ax.set_ylabel('Energy (eV)')
        ax.set_xlabel(f'Structure {structure}')
        ax.set_title(f'Energy Levels of GaInAsP/InP')
        ax.set_xticks([])  # Remove x-axis ticks
        ax.legend()

        # Annotate the energy gap between the valence band of the epilayer and the valence band of the substrate
        ax.text(0.51, (InP_Ev + GaInAsP_Ev) / 2, 
                f'Energy Gap: {(GaInAsP_Ev - InP_Ev):.2f} eV', 
                ha='center', va='center', fontsize=8, color='black')
        ax.annotate(f'',
                    xy=(0.4, InP_Ev) , xycoords='data',
                    xytext=(0.4, GaInAsP_Ev), textcoords='data',
                    arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"),
                    fontsize=8, ha='center', va='top')

        # Annotate the valence band and the conduction band
        ax.text(-0.2, InP_Ec, 'Ec', ha='right', va='center',
                fontsize=8, color='black')
        ax.text(-0.2, InP_Ev, 'Ev', ha='right', va='bottom',
                fontsize=8, color='black')
        ax.text(0.6, GaInAsP_Ec, 'Ec', ha='left', va='center',
                fontsize=8, color='black')
        ax.text(0.6, GaInAsP_Ev, 'Ev', ha='left', va='center',
                fontsize=8, color='black')
        ax.text(1.0, GaInAs_Ec, 'Ec', ha='left', va='center',
                fontsize=8, color='black')
        ax.text(1.0, GaInAs_Ev, 'Ev', ha='left', va='center',
                fontsize=8, color='black')

        plt.show()  # Display each subplot

def plot_GaInAs_InP(ax):
    """
    Plot the energy levels for GaInAs/InP structures.

    Parameters:
        ax (matplotlib.axes.Axes): The subplot axes.
    """
    # Get the energy levels of InP
    InP_Ec = krijin_table1.loc["InP", 'Ec']
    InP_Ev = krijin_table1.loc["InP", 'Ev']

    # Get T_ABC values for all structures for each material parameter
    structures = ingaas_table.columns[0:].tolist()
    GaInAs_Ec_list = process_data_ternary('Ec')
    GaInAs_Ev_list = process_data_ternary('Ev')

    # Iterate over structures
    for i, structure in enumerate(structures):
        GaInAs_Ec = GaInAs_Ec_list[i]
        GaInAs_Ev = GaInAs_Ev_list[i]

        # Create a new subplot for each structure
        fig, ax = plt.subplots()

        ax.bar([0], [InP_Ec - InP_Ev], bottom=InP_Ev, width=0.4, color='blue', label='InP')
        ax.bar([0.4], [GaInAs_Ec - GaInAs_Ev], bottom=GaInAs_Ev, width=0.4, color='red', label='GaInAs')
        ax.set_ylabel('Energy (eV)')
        ax.set_xlabel(f'Structure {structure}')
        ax.set_title(f'Energy Levels of GaInAs/InP')
        ax.set_xticks([])  # Remove x-axis ticks
        ax.legend()

        # Annotate the energy gap between the valence band of the epilayer and the valence band of the substrate
        ax.text(0.51, (InP_Ev + GaInAs_Ev) / 2, 
                f'Energy Gap: {(GaInAs_Ev - InP_Ev):.2f} eV', 
                ha='center', va='center', fontsize=8, color='black')
        ax.annotate(f'',
                    xy=(0.4, InP_Ev) , xycoords='data',
                    xytext=(0.4, GaInAs_Ev), textcoords='data',
                    arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"),
                    fontsize=8, ha='center', va='top')

        # Annotate the valence band and the conduction band
        ax.text(-0.2, InP_Ec, 'Ec', ha='right', va='center',
                fontsize=8, color='black')
        ax.text(-0.2, InP_Ev, 'Ev', ha='right', va='bottom',
                fontsize=8, color='black')
        ax.text(0.6, GaInAs_Ec, 'Ec', ha='left', va='center',
                fontsize=8, color='black')
        ax.text(0.6, GaInAs_Ev, 'Ev', ha='left', va='center',
                fontsize=8, color='black')

        plt.show()  # Display each subplot


def plot_quaternary(ax):
    """
    Plot the energy levels for quaternary alloys.

    Parameters:
        ax (matplotlib.axes.Axes): The subplot axes.
    """
    GaInAsP_Ec = process_data_quaternary('Ec')
    GaInAsP_Ev = process_data_quaternary('Ev')

     # Get T_ABC values for all structures for each material parameter
    structures = ingaas_table.columns[0:].tolist()
    GaInAs_Ec_list = process_data_ternary('Ec')
    GaInAs_Ev_list = process_data_ternary('Ev')
    GaInAs_Eg_list = process_data_ternary('Eg(Γ)')

    # Iterate over structures
    for i, structure in enumerate(structures):
        GaInAs_Ec = GaInAs_Ec_list[i]
        GaInAs_Ev = GaInAs_Ev_list[i]
        GaInAs_Eg = GaInAs_Eg_list[i]

        # Create a new subplot for each structure
        fig, ax = plt.subplots()

        # Iterate over structures and plot the values
        ax.hlines(y=GaInAsP_Ev, xmin=0, xmax=0.4, color='blue', label='Valence band')
        ax.hlines(y=GaInAsP_Ec, xmin=0, xmax=0.4, color='red', label='Conduction band')
        ax.vlines(x=0.4, ymin=GaInAsP_Ev, ymax=GaInAs_Ev, color='blue')
        ax.vlines(x=0.4, ymin=GaInAsP_Ec, ymax=GaInAs_Ec, color='red')
        ax.hlines(y=GaInAs_Ev, xmin=0.4, xmax=0.8, color='blue')
        ax.hlines(y=GaInAs_Ec, xmin=0.4, xmax=0.8, color='red')
        ax.vlines(x=0.8, ymin=GaInAsP_Ev, ymax=GaInAs_Ev, color='blue')
        ax.vlines(x=0.8, ymin=GaInAsP_Ec, ymax=GaInAs_Ec, color='red')
        ax.hlines(y=GaInAsP_Ev, xmin=0.8, xmax=1.2, color='blue')
        ax.hlines(y=GaInAsP_Ec, xmin=0.8, xmax=1.2, color='red')

        ax.set_ylabel('Energy (eV)')
        ax.set_xlabel(f'Structure {structure}')
        ax.set_title(f'Energy Levels for GaInAsP/GaInAs')
        ax.set_xticks([])  # Remove x-axis ticks
        ax.legend()

        ax.text(0.2, (GaInAsP_Ev + GaInAsP_Ec) / 2, 
                'GaInAsP', 
                ha='center', va='center', fontsize=8, color='black')
        ax.text(0.5, (GaInAs_Ev + GaInAs_Ec) / 2, 
                'GaInAs', 
                ha='center', va='center', fontsize=8, color='black')
        ax.text(1.0, (GaInAsP_Ev + GaInAsP_Ec) / 2, 
                'GaInAsP', 
                ha='center', va='center', fontsize=8, color='black')
        ax.text(0.7, (GaInAs_Ev + GaInAs_Ec) / 2, 
                f'Energy Gap: {GaInAs_Eg:.2f} eV', 
                ha='center', va='center', fontsize=8, color='black')
        ax.annotate(f'',
                    xy=(0.6, GaInAs_Ev), xycoords='data',
                    xytext=(0.6, GaInAs_Ec), textcoords='data',
                    arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"),
                    fontsize=8, ha='center', va='top')
        
        plt.show()  # Display each subplot


def main():

    # Get calculation type from user
    calculation_type = input("Enter the type of calculation (binary, ternary, quaternary): ")

    if calculation_type == 'binary':
        # Ask the user for the material combination (epilayer/substrate)
        while True:
            epilayer_material = input("Enter the epilayer material: ")
            substrate_material = input("Enter the substrate material: ")

            # Check if the entered materials are in the table
            if epilayer_material in krijin_table1.index and substrate_material in krijin_table1.index:
                break
            else:
                print("Invalid materials. Please enter valid materials.")
        
        # Process 001 Orientation
        fig, ax = process_data_binary('001', epilayer_material, substrate_material)

        # Process 111 Orientation
        fig, ax = process_data_binary('111', epilayer_material, substrate_material)

    elif calculation_type == 'ternary':
        fig, ax = plt.subplots()
        plot_GaInAs_InP(ax)
        plt.show()

    elif calculation_type == 'quaternary':
        fig, ax = plt.subplots()
        plot_GaInAsP_InP(ax)
        plt.show()

        fig, ax = plt.subplots()
        plot_quaternary(ax)
        plt.show()

    else:
        print("Invalid calculation type. Please enter 'binary', 'ternary', or 'quaternary'.")

if __name__ == "__main__":
    main()