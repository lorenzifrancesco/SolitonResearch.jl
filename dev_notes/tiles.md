##### Filling grid may be suicidal for armageddon-like simulations
simulations on the cloud confirm that it is.

##### 28/06/2023
##### We implement a serial version of the tiling procedure in order to not fill the GPU mem
- [x] implement the serial version of the tiling procedure
- [x] run it on the cloud -> tiles are all around $50\%$ transmission.
- [x] check $t_f$ -> it was all a matter of transforming the initial state

##### 3/07/2023
- [x] check tile by tile the "strange collapse" of the  3D GPE
- [x] design the values of velocity to find the critical barrier strength.
- [x] launch lines simulations.
How can we interpret the results of the 3D GPE and high nonlinearity?

##### Implement dynamic time_steps assignment for better quality tiles and lines
- [x] question: how can i estimate the right time_steps? I should have an almost constant $dt$ for having a good description of the phyisics!
- [x] rewrite the normal time axis scheduling such that $dt$ is constant -> eventually warn if time_steps are too many.
  - [x] briefly study the maximum $dt$ we can set for 3D and 1D simulations ->
     - 1D -> $dt=0.15$
     - 3D -> $dt=0.1$??

#### 3D tiles for $\gamma=0.65$ show some breaking at high velocity and barrier. 
What is the physical meaning?
- [ ] try to use a very precise $dt$ for that.

```
sd = load_parameters_alt(gamma_param=0.65)
       prepare_for_collision(sd, 0.65)
       g3 = sdc["G3"]

```