##### 3/07/2023
We need to obtain $\sigma^2$ for a single state (**ground state** or **collision center**) for 
- NPSE
  - [x] just run ```sigma2``` on the state
- NPSE_plus
  - [x] just solve the nonlinear problem ```sigma2``` on the state
- 3D GPE
  - [x] compute an estimate via least squares (variance estimator).
  the estimate routine is extremely slow. Vectorize it in GPU.
> the 3D case is difficult: the estimate seems not to be accurate due to the "sprayed values" at low axial probability.
- [x]  fix the 3D method: why is it at zero?

#### Comparison script?
we have a script for various values of the barrier. 
- [x] write a comparison line plotter