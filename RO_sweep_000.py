import numpy as np
from scipy.stats import skew, kurtosis
import sys

##################################################################################
### 1. RO solver                                                               ###
### Please manually uncomment or comment out white noise and red noise options ###
### inside the module. Default is white noise.                                 ###
##################################################################################

##################################################################################
### 1. RO solver                                                               ###
### Please manually uncomment or comment out white noise and red noise options ###
### inside the module. Default is white noise.                                 ###
##################################################################################

def simulate_sde_system(R0, Ra, F1, eps, F2, bT, cT, dT, bh, mT, sigmaT, B, mH, sigmaH, dt, T_total, T_initial, h_initial, xiT_initial, xih_initial, random_num, save_interval):

    ### Initializations ###
    num_steps = int(T_total / dt)
    t = np.arange(0, T_total, dt)
    T = np.zeros(num_steps)
    xiT = np.zeros(num_steps)
    h = np.zeros(num_steps)
    xiH = np.zeros(num_steps)

    T[0] = T_initial
    xiT[0] = xiT_initial
    h[0] = h_initial
    xiH[0] = xih_initial
    omega_a = 2*np.pi/12 

    np.random.seed(random_num)
    noiseT = np.random.normal(0, sigmaT * np.sqrt(dt), size=num_steps)
    np.random.seed(10000+random_num)
    noiseH = np.random.normal(0, sigmaH * np.sqrt(dt), size=num_steps)

    for i in range(1, num_steps):
        ### Introduce seasonality in R ###
        R = R0 + Ra * np.sin(omega_a * t[i])

        ### Introduce multiplicative noise where T > 0 ###
        if T[i-1] >= 0:
            my_B = B
        else:
            my_B = 0

        ### White noise for T ###
        xiT[i-1] = noiseT[i-1] 
        ### Red noise for T ###
        # xiT[i] = xiT[i-1] + dt * (-mT * xiT[i-1] + np.sqrt(2*mT) * noiseT[i-1]) 
        ### Update T ###
        T[i] = T[i-1] + dt * (R * T[i-1] + F1 * h[i-1]
                      + bT * T[i-1]**2   - cT * T[i-1] ** 3 
                      + dT * T[i-1] * h[i-1]) + xiT[i-1] * (1 + my_B * T[i-1]) 



        ### White noise for h ###
        xiH[i-1] = noiseH[i-1]
		### Red noise for h ###
        # xiH[i] = xiH[i-1] + dt * (-mH * xiH[i-1] + np.sqrt(2*mH) * noiseH[i-1])
        ### Update h ###
        h[i] = h[i-1] + dt * (-eps * h[i-1] - F2 * T[i-1] 
                      - bh * T[i-1]**2) + xiH[i-1] 
        

    return T[::save_interval], h[::save_interval]

##################################################################################
### 1. End of RO solver                                                        ###
##################################################################################




##################################################################################
### 2. Define parameters for sweep simulations                                 ### 
##################################################################################

### Linear parameters ###
R0_values = np.array([0.03, 0.05, 0.07]) 
Ra_values = np.array([0.16])
F1_values = np.array([0.14])
eps_values = np.array([0.13]) 
F2_values = np.array([0.16])

### Deterministic nonlinear parameters ###
bT_values = np.array([0.02])
cT_values = np.array([0.001])
dT_values = np.array([0.0])
bh_values = np.array([0.03])

### Noise parameters ###
sigmaT_values = np.array([0.20])
sigmaH_values = np.array([0.20])
B_values = np.array([0.4])

### Red noise parameters ###
mT_values = np.array([0.0])
mH_values = np.array([0.0])

##################################################################################
### 2. End of parameter definitions                                            ### 
##################################################################################



##################################################################################
### 3. Other simulation setups                                                 ### 
### The number "000" in the simulation code name RO_sweep_000.py refers to     ###
### random_num (the random seed number).                                       ###
##################################################################################

random_num = 0        # Random seed or reference number
T_total = 1920        # Total simulation length [months]
dt = 0.01             # Simulation time step [months]
T_initial = 1.0       # Initial condition for T
h_initial = 0.0       # Initial condition for h
xiT_initial = 0.0     # Initial value for T-noise
xih_initial = 0.0     # Initial value for h-noise
skip = 120            # Number of initial months to skip in statistical analysis
data_interval = 1.0   # Data saving interval [months]
save_interval = int(data_interval / dt)  # Save data every 'save_interval' steps

##################################################################################
### 3. End of other simulation setups                                          ### 
##################################################################################




##################################################################################
### 4. Perform parameter sweep simulations                                     ### 
### simul_{random_num}_{simul_num}.txt saves the time series of T and h,       ###
### while meta_{random_num}_{simul_num}.txt saves parameter information        ###
### and data statistics.                                                       ###
### simul_num (simulation reference number) is updated sequentially.           ###
##################################################################################

simul_num = 0      
for R0 in R0_values:
    for Ra in Ra_values:
        for eps in eps_values:
            for F1 in F1_values:
                for F2 in F2_values:
                    for bT in bT_values:
                        for cT in cT_values:
                            for dT in dT_values:
                                for bh in bh_values:
                                    for mT in mT_values:
                                        for sigmaT in sigmaT_values:
                                            for B in B_values:
                                                for mH in mH_values:
                                                    for sigmaH in sigmaH_values:
                                                        print(f"""random_num: {random_num}, simul_num: {simul_num}: 
                                                        R0: {R0}, Ra: {Ra}, F1: {F1}, F2: {F2}, eps: {eps}, 
                                                        bT: {bT}, cT: {cT}, dT: {dT}, bh: {bh}, 
                                                        mT: {mT}, sigmaT: {sigmaT}, B: {B}, mH: {mH}, sigmaH: {sigmaH}""")


														### Solve RO equations ###
                                                        T, h = simulate_sde_system(R0, Ra, F1, eps, F2, 
                                                                                   bT, cT, dT, bh, 
                                                                                   mT, sigmaT, B, mH, sigmaH, 
                                                                                   dt, T_total, T_initial, h_initial, 
                                                                                   xiT_initial, xih_initial, 
                                                                                   random_num, save_interval)

                                                        ### Calculate statistics ###
                                                        mean_T = np.mean(T[skip:])
                                                        std_T = np.std(T[skip:], ddof=1)
                                                        skew_T = skew(T[skip:])
                                                        kurtosis_T = kurtosis(T[skip:], fisher=False)

                                                        mean_h = np.mean(h[skip:])
                                                        std_h = np.std(h[skip:], ddof=1)
                                                        skew_h = skew(h[skip:])
                                                        kurtosis_h = kurtosis(h[skip:], fisher=False)

														### Save time series ###
                                                        data = np.array([T[skip:], h[skip:]]).T
                                                        np.savetxt(f"simul_{random_num:03d}_{simul_num:05d}.txt", (data), fmt='%.3f %.3f')

														### Save meta data ###
                                                        meta_data = np.array([mean_T, std_T, skew_T, kurtosis_T,
                                                                              mean_h, std_h, skew_h, kurtosis_h,
                                                                              random_num, simul_num,
                                                                              T_total, dt, data_interval,
                                                                              T_initial, h_initial, 
                                                                              xiT_initial, xih_initial, 
                                                                              R0, Ra, F1, eps, F2, 
                                                                              bT, cT, dT, bh, 
                                                                              mT, sigmaT, B, mH, sigmaH]).T


                                                        headers = ["mean(T)", "std(T)", "skew(T)", "kurt(T)", 
                                                                   "mean(h)", "std(h)", "skew(h)", "kurt(h)",
                                                                   "random_num", "simul_num",
                                                                   "T_total [months]", "dt [months]", "save_interval [months]",
                                                                   "T_initial", "h_initial", 
                                                                   "xiT_initial", "xih_initial",
                                                                   "R0", "Ra", "F1", "eps", "F2", 
                                                                   "bT", "cT", "dT", "bh", 
                                                                   "mT", "sigmaT", "B", "mH", "sigmaH"]

                                                        formatted_lines = ["{:>10.4f} # {}".format(data, header) 
                                                                           for data, header in zip(meta_data, headers)]

                                                        with open(f"meta_{random_num:03d}_{simul_num:05d}.txt", "w") as file:
                                                            file.write("\n".join(formatted_lines))
                                                    simul_num += 1

##################################################################################
### 4. End of parameter sweep simulations                                     ### 
##################################################################################



