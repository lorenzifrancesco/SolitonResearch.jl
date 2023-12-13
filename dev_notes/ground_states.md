# Ground state issues
- 1D ground state convergence fails
 
1D-GPE with $N=1024, \, a_{tol}=10^{-4}, \, dt = 0.01 $ 
---> no advantage in decreasing tolerance or $dt$.
 Error on $\mu$ of around $1\%$. The only sensible improvement is in increasing $N$.

- 3D ground state convergence is bound to the initial state, too much

- NPSE_plus smash the RAM

And we don't know why. Using GC.gc() helps but cost too much time.


##### Tiling 
Almost ok, but needs to be tested before the armageddon-like simulations.

#### TODO
- [x] 1D ground state convergence (rinuncio)
- [ ] Fix automatic ground state (wish)
- [ ] Fix CrankNicholson ground state (wish)



#### DO WE STILL HAVE CORRECTNESS PROBLEMS IN THE TRANSFORMS? 
 SEEMS SO, SINCE INCREMENTING THE STEPS IN THE SOLITON SIMULATION LEADS TO A DIFFERENT RESULT. _yes, but how do they SCALE?_
 - Parseval is well tested now.
- Maybe the transform domain is somehow corrupted and not taking into account the right nonlinearity?

---> working in kspace, we are "inverted"

### 27/06/2023
### Debug of strange behaviour of ground states
- [x] TEST: CrankNicholson -> seems to be without kinetic term: check initialization of tridiagonal matrices
- [x] - TEST: CrankNicholson initialized with soliton solution -> broken (seems without kinetic term)
- [x] TEST: SplitStep initialized with soliton solution -> After relatively few iterations, it converges to the wrong function (slightly below)  
- [x] BackwardEuler is shifting the ground state (not correct)

- [ ] Quick FIX CrankNicholson (nota available)

- [ ] Check correctness of the tranforms


**NB: we still do not have 3D CrankNicholson**
##### Curiosity
- [ ] obtain the Thomas-Fermi function simulating without kinetics

#### Why SplitStep is broken? Transforms??
We should check the kspace, the xspace should be ok (the transforms are inverse of each other)
- [x] check using $dk = 2\pi/L$. It is used like that! -> using it in the transforms breaks everything.
- [ ] justify theoretically check the correct version for $dk$ using sampling arguments

**NB: in our code the multiplicative constant to have the correct sampling of the transform are written in a counterintuitive way (DX, DK)**.

- in measures() the $dk$ is not corresponding to $2/pi$
- [ ] check also chempot functions

#### we are using a naive kvec!! 
- [x] try to replace the correct version, the tests survives
- [x] try to see if something changes in the soliton simulation: NOTHING CHANGES

There are two different $dk$, one used for norm in kspace - **modified** ($dk=\frac{2\pi}{dx \,N }$), and one which is the official one of the discretization **official** ($dk = \frac{2\pi}{L}$). 
WHY?

- the one used in kvecs (official) affects also the phyisics: it is contained in k^2. 
- In the FFT the only effect of the choice of $dk$ is the prefactor chosen for the normalization of the transform. In such a factor the "$dk$" is the modified

- [x] can we modify the kvec in the FFT to have the correct $dk$? without breaking down everything?

##### find the scaling of the error
Error scales strangely like $\sim 1/N^2$. This is expected for the SplitStep ??
- [ ] check in the literature: is this scaling expected? (curiosity)

_NB: test still fails with 1D GPE GS_
##### 28/06/2023
#### MACRO problem: convergence is too much dependent on the initial condition.

##### TODO - Fix of the 3D covergence
The 3D result stays too close to the initial state in the axial profile. Maybe it adjusts in the radial dimension, but how?
> we work in a seqential way using as initial conditions the ground states of the previous simulations.
- [x] focus the problem in gamma = 0.6 case -> best suited script: soliton_gs_comparison.jl

![](2023-06-28-13-52-22.png)
![](2023-06-28-13-52-50.png)
- [x] try not to optimize the initial condition
  - [x] try a flat one
  ![](2023-06-28-13-53-30.png)
  - [x] try a sharp one
  ![](2023-06-28-13-53-00.png)
  Essentially, the convergence routine for the 3D case is BROKEN (try to see the initial states together? nah).
  **NB results with $\text{abstol} = 10^{-3}$**.
- [x] try to use different tolerances. Maybe we get to see the convergence. Test with $\text{abstol} = 10^{-5}$ -> nothing changes in _sharp_ case.   
Also changing N does not help.
- [x] What about changing N_transverse? NOTHING (flat case)
![](2023-06-28-14-31-38.png)
we double 64 -> 128. Very slow iteration speed. 
- [x] visualize the transverse distribution (may give a hint).
-> prepared visualization functions. The visualization of the wrong results does not suggest strange behaviours.


##### Make the method converge
- [x] check for a quick fix -> QUICK FIX: we were plotting the initial condition.

##### Other strange issue
When we copy the initial state to the 3D routine, then we 
This could suggest us the reason behind the difference between the analytical and numerical 1D-GPE.

It may be the same problem or not.

- [ ] check for a quick fix.

## This is a problem even for the transverse solution of the 3D GPE. Imprecision of a few percents...
- [ ] debug it! 

MWE