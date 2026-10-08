# Beware that this does not use the common cpp header 
# file, so if there is a change in it, this code has to be 
# updated manually.

import numpy as np
import matplotlib.pyplot as plt

kT_center = np.array([0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65]);
energies = np.array([3.0, 3.2, 3.5, 3.9, 4.5, 7.7, 9.2, 11.5, 14.5, 19.6, 27.0]);

rmaxdef = 80.0#72.0

def rhofitmax(kT, energy, syst=0):
    ktparam = (kT-0.2) / 0.05
    #ktparam = np.sqrt(ktparam) # tried but this goes down too slow, even if multiplied not by 5.0 but by 10.0
    systpm = 0.0
    if syst == 1:
        systpm = -1.0
    elif syst == 2:
        systpm = 1.0
    return (rmaxdef - ktparam * 5.0) * (1.0 + systpm * 0.25) * np.sqrt(energy/11.5) # only pm20% instead of the previous 60%

plt.rcParams['axes.titlesize'] = 20   # Title font size
plt.rcParams['axes.labelsize'] = 20   # Axis labels font size
# Also set larger fontsize for legends and ticks:
plt.rcParams['legend.fontsize'] = 16   # Legend font size
plt.rcParams['xtick.labelsize'] = 16   # X-axis tick labels font size
plt.rcParams['ytick.labelsize'] = 16   # Y-axis tick labels font size
# Do not leave white space on the sides of the plot:
plt.rcParams['figure.subplot.left'] = 0.05
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.top'] = 0.9
plt.rcParams['figure.subplot.bottom'] = 0.15


#plt.figure(figsize=(8,6))
#plt.title(r'$\rho_{fit}^{max}$ vs $k_T$ for different energies')
#plt.xlabel(r'$k_T$ (GeV)')
#plt.ylabel(r'$\rho_{fit}^{max}$')
#for energy in energies:
#    plt.plot(kT_center, rhofitmax(kT_center, energy), 'o-', label=f'{energy} GeV')
#plt.legend()
#fig1 = plt.gcf()
#fig1.savefig("figs/rhofitmax_vs_kt_snn.png")
#plt.show()

# Plot all three syst variations side-by-side
plt.figure(figsize=(21,6))
#plt.suptitle(r'$\rho_{fit}^{max}$ vs $k_T$ for different energies and systematic variations')
for syst in range(3):
    plt.subplot(1, 3, syst+1)
    if syst == 0:
        plt.title('Default')
    elif syst == 1:
        plt.title(r'Smaller $\rho_{fit}^{max}$')
    else:
        plt.title(r'Larger $\rho_{fit}^{max}$')
    plt.xlabel(r'$K_T$ (GeV)')
    plt.ylabel(r'$\rho_{fit}^{max}$')
    for energy in energies:
        plt.plot(kT_center, rhofitmax(kT_center, energy, syst), 'o-', label=f'{energy} GeV')
plt.legend()
fig2 = plt.gcf()
fig2.savefig("figs/rhofitmax_vs_kt_snn_all_syst.png")
#plt.show() # do not show if called from command line, to avoid blocking the script from finishing and saving the figure