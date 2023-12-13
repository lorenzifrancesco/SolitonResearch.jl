#### COLLAPSES
- The collapse gamma found for NPSE is found to be around 0.751, strange enough
This behavious depends strongly on the spatial discretization.
N=1024 -> 0.73qualcosa
N=2048 -> 0.751
N=4096 -> 0.74533

- NPSE_plus seems to have a higher threshold for collapse.
N=1024 -> gamma_c>0.75 
> (are we sure we are even collapsing at all in NPSE_plus?) 


##### TODO - Understand the strange collapse behaviour
- [x] try different N scales, and different initial conditions -> **nope**, only slight changes
- [x] try different number of steps: by increasing the simulation time, we may be able to see the collapse -> **nope**
##### what we do dyamical simulations
> The meaning of the parameter $\gamma_c$ is that for $\gamma > \gamma_c$ there cannot be a stable solitonic solution. 

Our method for finding solitonic solutions may not be probing their long-term stability. 
One may augment the precision in order to keep the simulation running for more time.

**Collapse script gives segmentation fault!!**
Big problem for prototyping

##### $\gamma_c$ is even bigger in the case of a dynamic simulation!!


##### 29/06/2023
##### we can also consider to have a comparison with the analytical chemical potential. (still not done!)
Intuitively we are ok with the chemical potential given. But we may check explicitly that the chemical potential is the same as the one given by the analytical solution.
- [x] we can check the analytical solution for the NPSE.
-> not easy, the implicit expression doesn't behave well

- [x] run some test at $\gamma = 2/3$ to see if we can reach collapse for some values

##### we may have problems with sigma2 method
- we do: when we change gamma we should change also sigma2 through init_sigma_2.
-> no, we have the correct sigma2 inside the simulations.
So how is it possible that we don't collapse in the simulations??

> issue found: ground state runs smoothly, but running sigma2 over the result panics a NPSE collapse. How is this possible? Maybe we cannot redefine effectively the sigma2 function? -> no, we just need to remember to transform back in xspace.

- we have a way at least to find the maximum value of the solution for a particular gamma:
    - [x] implement it
    - [x] check it against the plots in the paper: correct!
 - [x] (We can also check better our ground states with respect to the paper)

#### How does the NPSE-collapse suggest the true collapse?
We need to have a solution from below: 
- [x] maybe investigating the 3D collapse first would give us a better view of what a true collapse is.
- Currently, our best method for finding the 3D collapse is to stop the simulation when the probability of a point is very high. Let's say $0.05$.
- [x] NPSE-collapse can be reached "sharping up" a solution.

##### Lets try to find the collapse in the 3D case
-> introduce the Gpe3DCollapse exception: it seems to work well

##### Go on with NPSE simulations
- [x] try playing with the accuracy -> increasing the accuracy we seem to pinpoint lesser and lesser $\gamma_c$.
- [x] do a super accurate $dt=0.001$ $abstol=10^{-8}$ simulation -> we finally have an accurate collapse measurement. Commited as **sim** -> result $\gamma_c = 0.6698$
Maybe the problem is just the discretization (MAYBE).

- [x] check the NPSE solutions very near $\gamma_c$. Are they correct with respect to the maximum value?

##### Now we can focus on the 3D case!!

### Alternative way (faster): start from the ground state and see if the solver moves from that. If it does, it is going to collapse.