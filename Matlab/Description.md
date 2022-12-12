# Code Description
This folder stores the finite-difference code used for creating the Dataset related to the Duffing oscillator. 

The script `DuffingOSC` is used for computing a linear and a Duffing Oscillator with a linearly implicit scheme.
The script `DuffingOscSeparated` computes the Duffing oscillator with the same linearly implicit scheme seen in the other script. In addition, it computes the Duffing oscillator by calculating separately the linear and nonlinear parts at each time sample. The next sample is than given by the sum of the linear and nonlinear components. Beware that the solutions obtained with the two approaches differ by an error of the order of $10^{-9}$. The error seems to be linked to machine accuracy, because it shows random variations.

# Dataset Description
The datasets are made of 25 different solutions of the oscillators, with initial displacements and velocities betwen $[-1,1]$, computed with parameters $\omega = 1$, $\gamma = 1$, at a sample rate of 1000 Hz, for 120.1 sec.

1. outDOFull stores the Duffing oscillator solution
2. outSHO stores the linear oscillator solution
3. outLin stores the linear part of the duffing oscillator
4. outNlin stores the nonlinear part of the duffing oscillator
4. timeVec stores the time vector in samples
5. timeVecSec stores the time vector in seconds; i.e. timeVec*k, with k being the sampling step