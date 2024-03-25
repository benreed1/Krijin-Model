import math
import pandas as pd
import matplotlib.pyplot as plt

# Read data from files
krijin_table1 = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Krijin Table 1.xlsx', index_col='Material')
krijin_table2 = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Krijin Table 2.xlsx', index_col='Material')
data = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\data.xlsx', index_col='Material')
Data_InGaAsSb = pd.read_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Data_InGaAsSb.xlsx', index_col='Material')

def calculate_bandgap(substrate, epilayer, orientation):
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

    return Ev_epilayer, Ec_epilayer, Ev_substrate, Ec_substrate, ε_parallel, ε_perpendicular

def process_data_binary(orientations, epilayer_material, substrate_material):

    # Get material parameters for the specified epilayer and substrate
    epilayer_data = krijin_table1.loc[epilayer_material].values
    substrate_data = krijin_table1.loc[substrate_material].values

    # Calculate bandgap energies for orientations
    Ev_epilayer, Ec_epilayer, Ev_substrate, Ec_substrate, ε_parallel, ε_perpendicular = calculate_bandgap(substrate_data, epilayer_data, orientation=orientations)

    # Compile data into spreadsheet
    data.loc[epilayer_material, 'Ev'] = Ev_epilayer
    data.loc[epilayer_material, 'Ec'] = Ec_epilayer
    data.loc[epilayer_material, 'ε_parallel'] = ε_parallel
    data.loc[epilayer_material, 'ε_perpendicular'] = ε_perpendicular
    data.loc[substrate_material, 'Ev'] = Ev_substrate
    data.loc[substrate_material, 'Ec'] = Ec_substrate
    data.to_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\data.xlsx')

def calculate_T_values_quaternary(GaAs, InAs, GaInAs, InSb, GaInSb, InAsSb, GaSb, GaAsSb, material_parameter):
    x = 0.14
    y = 0.3
    if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0':
        T_ABC_values_quaternary = x * GaAs + (1 - x) * InAs + x * (x - 1) * GaInAs
        T_ABD_values_quaternary = x * InSb + (1 - x) * GaSb + x * (x - 1) * GaInSb
        T_ACD_values_quaternary = y * GaAs + (1 - y) * InSb + y * (y - 1) * InAsSb
        T_BCD_values_quaternary = y * InAs + (1 - y) * GaSb + y * (y - 1) * GaAsSb
    elif material_parameter == 'Ev' or material_parameter == 'Ec' or material_parameter == 'ε_parallel' or material_parameter == 'ε_perpendicular':
        T_ABC_values_quaternary = x * GaAs + (1 - x) * InAs 
        T_ABD_values_quaternary = x * InSb + (1 - x) * GaSb 
        T_ACD_values_quaternary = y * GaAs + (1 - y) * InSb 
        T_BCD_values_quaternary = y * InAs + (1 - y) * GaSb


    return T_ABC_values_quaternary, T_ABD_values_quaternary, T_ACD_values_quaternary, T_BCD_values_quaternary

def calculate_Q_values(T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values):
    # Calculate Q(x, y) using Equation (10)
    x = 0.14
    y = 0.3
    Q_values = (
        (x * (1 - x) * (y * T_ABC_values + (1 - y) * T_ABD_values)) / (x * (1 - x) + y * (1 - y)) +
        (y * (1 - y) * (x * T_ACD_values + (1 - x) * T_BCD_values)) / (x * (1 - x) + y * (1 - y))
    )

    return Q_values

def process_data_quaternary(material_parameter):
    # Get material parameters
    GaAs = krijin_table1.loc["GaAs", material_parameter]
    InAs = krijin_table1.loc["InAs", material_parameter]
    InSb = krijin_table1.loc["InSb", material_parameter]
    GaSb = krijin_table1.loc["GaSb", material_parameter]
    if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0':
        GaInAs = krijin_table2.loc["GaxIn1-xAs", f'C({material_parameter})']
        GaAsSb = krijin_table2.loc["GaAsxSb1-x", f'C({material_parameter})']        
        InAsSb  = krijin_table2.loc["InAsxSb1-x", f'C({material_parameter})'] 
        GaInSb = krijin_table2.loc["GaxIn1-xSb", f'C({material_parameter})'] 
        T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values = calculate_T_values_quaternary(GaAs, InAs, GaInAs, InSb, GaInSb, InAsSb, GaSb, GaAsSb, material_parameter)
    if material_parameter == 'Ev' or material_parameter == 'Ec' or material_parameter == 'ε_parallel' or material_parameter == 'ε_perpendicular' :
        GaInAs = GaInSb = InAsSb = GaAsSb = None
        T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values = calculate_T_values_quaternary(GaAs, InAs, GaInAs, InSb, GaInSb, InAsSb, GaSb, GaAsSb, material_parameter)
    # Calculate Q values for quaternary alloys
    Q_xy_values = calculate_Q_values(T_ABC_values, T_ABD_values, T_ACD_values, T_BCD_values)
    # Display the results
    if material_parameter == 'Eg(Γ)' or material_parameter == 'Δ0' or material_parameter == 'Ev' or material_parameter == 'Ec':
        print(f'Q(x,y): {material_parameter}: {Q_xy_values} eV ')
    elif material_parameter == 'ε_parallel' or material_parameter == 'ε_perpendicular':
        print(f'Q(x,y): {material_parameter}: {Q_xy_values}')

    # Compile data into spreadsheet
    if material_parameter == 'Eg(Γ)':
        Data_InGaAsSb.loc['InGaAsSb', 'Eg(Γ)'] = Q_xy_values
    elif material_parameter == 'Δ0':
        Data_InGaAsSb.loc['InGaAsSb', 'Δ0'] = Q_xy_values
    elif material_parameter == 'Ev':
        Data_InGaAsSb.loc['InGaAsSb', 'Ev'] = Q_xy_values
    elif material_parameter == 'Ec':
        Data_InGaAsSb.loc['InGaAsSb', 'Ec'] = Q_xy_values
    elif material_parameter == 'Ec':
        Data_InGaAsSb.loc['InGaAsSb', 'ε_parallel'] = Q_xy_values
    elif material_parameter == 'Ec':
        Data_InGaAsSb.loc['InGaAsSb', 'ε_perpendicular'] = Q_xy_values
    Data_InGaAsSb.to_excel(r'C:\Users\beree\OneDrive\Documents\ENGR301\Krijin Paper\Data_InGaAsSb.xlsx')

    return Q_xy_values
    

def main():

    # Get calculation type from user
    calculation_type = input("Enter the type of calculation (binary, quaternary): ")

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

        process_data_binary('001', epilayer_material, substrate_material)

    elif calculation_type == 'quaternary':
        #Run calculations to import data into spreadsheet
        process_data_quaternary('Eg(Γ)')
        process_data_quaternary('Δ0')
        process_data_quaternary('Ev')
        process_data_quaternary('Ec')
        process_data_quaternary('ε_parallel')
        process_data_quaternary('ε_perpendicular')

    else:
        print("Invalid calculation type. Please enter 'binary', or 'quaternary'.")

if __name__ == "__main__":
    main()