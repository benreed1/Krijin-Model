import matplotlib.pyplot as plt

WW = [8.5, 5, 3.6, 2.8, 2.3, 2]
x = [0.48, 0.39, 0.32, 0.25, 0.17, 0.15]
HHmin = [-0.48103464, -0.47558744, -0.46954987, -0.46144521, -0.44899266, -0.44557311]
HHmax = [-0.42817994, -0.39146782, -0.36531519, -0.34075079, -0.31393487, -0.30766182]
Emin = [0.31548451, 0.28501109, 0.26089222, 0.23656524, 0.20881289, 0.2016163]
Emax = [0.66338614, 0.66883335, 0.67487094, 0.68297562, 0.69542824, 0.69884781]
Ee = [0.37463068, 0.39459382, 0.41149551, 0.42341543, 0.42759945, 0.44287934]
Eh = [-0.45428465, -0.44559327, -0.44101724, -0.43544518, -0.42529846, -0.42572863]

for i in range(len(WW)):
    fig, ax = plt.subplots()

    ax.hlines(y=HHmin[i], xmin=-WW[i], xmax=-WW[i]/2, color='blue')
    ax.hlines(y=Emax[i], xmin=-WW[i], xmax=-WW[i]/2, color='red')
    ax.vlines(x=-WW[i]/2, ymin=HHmin[i], ymax=HHmax[i], color='blue')
    ax.vlines(x=-WW[i]/2, ymin=Emax[i], ymax=Emin[i], color='red')
    ax.hlines(y=HHmax[i], xmin=-WW[i]/2, xmax=WW[i]/2, color='blue')
    ax.hlines(y=Emin[i], xmin=-WW[i]/2, xmax=WW[i]/2, color='red')
    ax.hlines(y=Ee[i], xmin=-WW[i]/2, xmax=WW[i]/2, color='green')
    ax.hlines(y=Eh[i], xmin=-WW[i]/2, xmax=WW[i]/2, color='green')
    ax.vlines(x=WW[i]/2, ymin=HHmin[i], ymax=HHmax[i], color='blue')
    ax.vlines(x=WW[i]/2, ymin=Emax[i], ymax=Emin[i], color='red')
    ax.hlines(y=HHmin[i], xmin=WW[i]/2, xmax=WW[i], color='blue')
    ax.hlines(y=Emax[i], xmin=WW[i]/2, xmax=WW[i], color='red')

    ax.set_ylabel('Energy (eV)')
    ax.set_xlabel(f'Well Width')
    ax.set_title(f'Energy Levels for GaInAsP/GaInAs with Well Width={WW[i]}')
    ax.set_xlabel('L(m^-9)')
    ax.legend()

    #ax.text(-1.25*WW[i], (HHmin[i] + Emax[i]) / 2, 
            #'GaInAsP', 
            #ha='center', va='center', fontsize=8, color='black')
    #ax.text(-WW[i]/4, (HHmax[i] + Emin[i]) / 2, 
            #'GaInAs', 
            #ha='center', va='center', fontsize=8, color='black')
    #ax.text(1.25*WW[i], (HHmin[i] + Emax[i]) / 2, 
            #'GaInAsP', 
            #ha='center', va='center', fontsize=8, color='black')
    ax.text(WW[i]*0.3, (HHmax[i] + Emin[i]) / 2, 
            f'Energy Gap: {Ee[i] - Eh[i]:.2f} eV', 
            ha='center', va='center', fontsize=8, color='black')
    ax.annotate(f'',
                xy=(0, Ee[i]), xycoords='data',
                xytext=(0, Eh[i]), textcoords='data',
                arrowprops=dict(arrowstyle="<->", connectionstyle="arc3"),
                fontsize=8, ha='center', va='top')
plt.show()