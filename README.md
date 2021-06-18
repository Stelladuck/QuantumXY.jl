# quantumxy

Julia scripts for a 1D quantum spin-1/2 xy chain with periodic condition.
It contains the code for the Hamiltonian, its magnetization and correlation.
The scripts are written for the purpose of assignment.

## Functions
```generate_Hamiltonian(N, J, h)``` : Function generating Hamiltonian using basis\
```generate_Hamiltonian_with_Pauli_matrices(N, J, h)``` : Function generating Hamiltonian using Pauli matrices\
```magnetization(state, basis)``` : Function computing magnetization\
```correlation(state, basis, r)``` : Function computing correlation\
$N$ is the system size. $J$ is the coupling strength. $h$ is the magnetic field. $r$ is the degree of neighboring sites.
If $r=1$, it considers the neareset neighbors.
