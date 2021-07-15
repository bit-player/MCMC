# MCMC

Five variants of Monte Carlo simulations of the Ising model of ferromagnetism, all written in JavaScript and meant to run in a browser window.

This code is associated with the article "Three Months in Monte Carlo" at [http://bit-player.org/2021/three-months-in-monte-carlo](http://bit-player.org/2021/three-months-in-monte-carlo).

- Program 1: Metro vs. Glauber. Comparison of the Metropolis algorithm and Glauber dynamics.

- Program 2: Mix and Match. The Metropolis and Glauber algorithms differ in two main ways: the order in which sites are visited in the lattice, and the rule applied to decide whether or not a selected spin will be flipped. Program 2 allows those components to be recombined in the four possible ways.

- Program 3: Visitation Variations. Explore the effects of eight choices for the visitation sequence in the Metropolis algorithm.

- Program 4: Boundaries. In a computer model we can create only a finite chunk of what might be an infinite plane. Here are eight ideas for what the model might do when you come to the edge of the world.

- Program 5: The MCMC Microscope. A slow-motiom, close-up view of how individual lattice sites evolve under various rules and protocols.

Also included here is a Julia file (meant to be opened and run in the Pluto.jl notebook server) with data and programs for reproducing the graphs in the bit-player.org article.
