# Optimization-of-Levi-function-using-PSO

This project aims to estimate the Levi function's minimum value using Particle Swarm Optimization (PSO). 

PSO is an optimization algorithm used to decipher the global minimum and maximum global value of a non-linear continuous function. It is based on the concept of swarm intelligence capable of solving mathematical problems at engineering. 
Advantage: it has very few parameters that we would have to adjust. 

In the MATLAB code (POS) Section I, 
I have considered a swarm of 100 particles. 
Number of particles = 100
iterations = 100
Range of x and y fixed [-10,10]

Randomly assign the values of those particles in the given region. 
Initial velocities of all the 100 particles in 'x' and 'y' direction= 0

Particles (matrix)
1st column: position in 'x' direction
2nd column: position in 'y' direction
3rd column: velocity in particle in 'x' direction
 4th column: velocity in particle in 'y' direction 

1st iteration
Calculate all the particles' fitness values and put them in the 5th column of the 'particles' matrix.
Sort the rows in correspondence with the fitness values that are in the 5th column.
Set Pbest (100 particles x 3 - 'x' position, 'y' position, and fitness value) and Gbest, where Pbest is the personal best of each particle, and Gbest is the group best for each iteration.
Initially, Pbest would be a copy of the initial positions itself. Gbest is the best-fit particle in the 1st iteration. (Put this Pbest column in the particle-matrix only)
Make the Gbest list to keep track of Gbest(s) with iterations.

Now, update the velocities and the positions using the update rule equations.  

The next velocity of a particle depends on 
The previous velocity
The difference between Pbest and the position earlier 
And the difference between Gbest's position until that iteration and the position earlier. In this way, all the particles in the swarm can share information on the best point achieved regardless of which particle has found it.  

Parameters to be chosen in the velocity update rule:
Alpha - inertia weight constant, range - [0,1.2]
If alpha lies between [1, 1.2], it implies its previous motion entirely influences its motion. So, it may keep on going in the same direction. However, if alpha lies between [0,1), it increases the diversification, hence finding the global minimum increase (but simulation time would increase).

c1 (individual cognition parameter) - it raises the importance of the particle's own previous experience. range - [0,2]
r1 (random parameter) should be within [0,1] to avoid premature convergence. 

c2 (social learning parameter): it raises the importance of global learning of the swarm.
r2 (random parameter) should be within [0,1] to avoid premature convergence.  

The next position of the particle depends on:
its previous position
its new/next velocity

After this, if I find that some of the particles' positions exceed this range [-10,10], I can set those positions to the boundary values. 
Evaluate the Levi function's fitness values, followed by the sorting of rows w.r.t the fitness values. Compute the Gbest for this iteration and compare it with the previous value of Gbest. If it comes out to be lesser than the previous one (since this is a minimization problem) then, replace the Gbest value otherwise don't. Fill the selected Gbest into the 'Gbest_list.' For updating the Pbest values, if the new Pbest is lesser than the previous one, then update it otherwise don't. 

This procedure is repeated for 100 iterations.  
In the end, the Gbest value is the answer to the minimization problem.


In the MATLAB code (POS) Section II,
An animated graph is created where the Levi function is plotted, and the Pbest values for the first iteration are marked. 

In the MATLAB code (POS) Section III, 
An animated graph is created that displays how the Gbest moves over the surface of the Levi function.

In the MATLAB code (POS) Section IV, 
a graph is plotted between the number of iterations (x-axis) and the fitness values of Gbest positions (y-axis).

In the MATLAB code (Pbest_Plot_Convergence),
An animated graph is created showing the convergence of Pbest values, considering all the iterations.


References: 
https://www.youtube.com/watch?v=nnkTSX5U_a4
https://www.intechopen.com/books/swarm-intelligence-recent-advances-new-perspectives-and-applications/particle-swarm-optimization-a-powerful-technique-for-solving-engineering-problems

