# Maxwell-Demon-in-open-quantum-systems
In this project I´ll explore a quantum implementation of a Maxwell’s demon using coupled quantum dots and open quantum system dynamics by means of density matrix and Lindblad formalisms. QuTiP-based simulations and codes, will be added in future iterations. As for now, I´m only uploading the thesis itself and a couple of calculators. Also, other demons will be studied as well


## Project status and future directions

This repository currently contains selected results and a first working implementation
developed for my Bachelor’s Thesis on a quantum Maxwell’s demon implemented with
coupled quantum dots and phonon/electronic baths. I recommend reading the whole .pdf for a full explanation on the topic.

The present code corresponds to a minimal and functional version of the model:
- Diagonalization of the interacting Hamiltonian (finding egenstates basis)
- Electronic and phononic transition matrices
- Liouvillian construction (finding superoperator)
- Steady–state populations
- Electronic, phononic and total heat currents
- Power and entropy flows
- Identification of thermodynamic regimes (engine, refrigerator and heat pump)
- Discussing the demon `ಠ‿ಠ´ψ

This is **not the final version of the project** (I´ll be working on it from time to time), but an initial base that will be significantly
extended in the following directions:

- Migration to **QuTiP** for a more systematic and scalable open quantum system treatment. Not only thermal demons are also desirable to study
- Exploration of additional system geometries and parameter regimes
- Study of information–thermodynamics relations and Maxwell-demon-like feedback
- Extension towards light–matter and cavity-coupled configurations

The goal is to develop a flexible and extensible framework to have fun with
non-equilibrium quantum thermodynamics in nanoscale systems.
