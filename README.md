README
-----------
README: RO Parameter Sweep Simulations (v0.0, May 12, 2025)
Author: Sooman Han (sooman.han@yale.edu)


DESCRIPTION
-----------
This project performs parameter sweep simulations of a Recharge Oscillator (RO)
model to study key ENSO dynamics, including amplitude asymmetry, duration
asymmetry, and regime classification (damped vs. self-sustained). The
Euler–Maruyama method is used to numerically integrate the stochastic
differential equations in the RO model, accounting for both additive and
multiplicative noise components.

The code is designed with minimalistic logic to efficiently explore a wide
parameter space, rather than relying on a single parameter set obtained through
conventional fitting methods.

Each simulation generates:
- A time series of temperature (T) and thermocline depth (h) anomalies
- Metadata including parameters used and summary statistics

**Please cite:**
Han, S., A. V. Fedorov and J. Vialard, 2025: *Realistic ENSO Dynamics Requires
a Damped Nonlinear Recharge Oscillator*. Submitted to *Journal of Climate*.

The code is also used in: 
Fangyu, F. J. Vialard, S. Han, Y. Planton, M. Lengaigne, G. Srinivas, S. Zhao, 
E. Guilyardi, W. Zhang, C. Ethé, R. Person, A. Voldoire, F.-F. Jin,
A. V. Fedorov and M. J. McPhaden, 2025: *ENSO cycles mostly after extreme 
El Niño events*. Submitted to *Nature*.

--------------------------------------------------------------------------------
FILES
-----

1. **RO_sweep_{random_num}.py**
   - Simulation scripts using a different random seed (`random_num`).

2. **RO_generator.sh**
   - Bash script to generate `RO_sweep_{random_num}.py` files from `RO_sweep_000.py`.

3. **simul_{random_num}_{simul_num}.txt**
   - Time series output of T and h.
   - Format: 2 columns (T, h)

4. **meta_{random_num}_{simul_num}.txt**
   - Contains parameter values and computed statistics (e.g., mean, std, skewness).

--------------------------------------------------------------------------------
HOW TO RUN
----------

1. The only required Python dependencies are `numpy` and `scipy`.

2. Set up your simulation by editing `RO_sweep_000.py`.  
   *Do not change the filename or the line `random_num = 0` inside this file.*

3. To generate multiple sweep scripts, open `RO_generator.sh` and set:

   ```bash
   start=1
   end=9  # for example

4. Then run the script to generate multiple Python simulation files:

   ./RO_generator.sh

   This will generate RO_sweep_001.py to RO_sweep_009.py,
   each with random_num = 1 to random_num = 9 internally.

5. Run each simulation script individually:

   ```bash
   python RO_sweep_000.py
   python RO_sweep_001.py
   ...
   python RO_sweep_009.py
