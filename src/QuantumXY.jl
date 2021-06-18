module QuantumXY

using IterTools
using LinearAlgebra

export generate_full_basis
export generate_Hamiltonian
export generate_Hamiltonian_with_Pauli_matrices
export magnetization
export correlation

# Function generating full basis
function generate_full_basis(N)
    return reshape(collect(IterTools.product([[1,0] for i=1:N]...)), 2^N)
end
# Example:
# generate_full_basis(3)
# 8-element Array{Tuple{Int64,Int64,Int64},1}:
#  (1, 1, 1)
#  (0, 1, 1)
#  (1, 0, 1)
#  (0, 0, 1)
#  (1, 1, 0)
#  (0, 1, 0)
#  (1, 0, 0)
#  (0, 0, 0)


# Function generating Hamiltonian using basis
function generate_Hamiltonian(basis::Vector{Tuple{Int64, Int64, Int64}}, J, h)
    M = length(basis) # Number of basis
    N = length(basis[1]) # Number of sites
    H = zeros(Float64, M,M) # Initialization of Hamiltonian
    # Update Hamiltonian
    for i in 1:M
        for j in 1:N
            if j < N
                # Update flipping terms
                swapped = collect(basis[i])
                if swapped[j] != swapped[j+1]
                    swapped[j], swapped[j+1] = swapped[j+1], swapped[j]
                    H[i, findfirst(x -> x == tuple(swapped...), basis)] += J
                elseif swapped[j] == swapped[j+1]
                    H[i, findfirst(x -> x == tuple(swapped...), basis)] = 0
                end
            elseif j == N
                # Periodic boundary condition
                swapped = collect(basis[i])
                if swapped[j] != swapped[1]
                    swapped[j],swapped[1] = swapped[1],swapped[j]
                    H[i, findfirst(x -> x == tuple(swapped...), basis)] += J
                elseif swapped[j] == swapped[1]
                    H[i, findfirst(x -> x == tuple(swapped...), basis)] = 0
                end
            end
        end
        # Update external field terms
        H[i, i] += h*(count(x -> x == 1, basis[i]) - count(x -> x == 0, basis[i]))
    end
    return -H
end
# Example:
# 8×8 Matrix{Float64}:
#  -3.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0
#  -0.0  -1.0  -1.0  -0.0  -1.0  -0.0  -0.0  -0.0
#  -0.0  -1.0  -1.0  -0.0  -1.0  -0.0  -0.0  -0.0
#  -0.0  -0.0  -0.0   1.0  -0.0  -1.0  -1.0  -0.0
#  -0.0  -1.0  -1.0  -0.0  -1.0  -0.0  -0.0  -0.0
#  -0.0  -0.0  -0.0  -1.0  -0.0   1.0  -1.0  -0.0
#  -0.0  -0.0  -0.0  -1.0  -0.0  -1.0   1.0  -0.0
#  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0   3.0


# Function generating Hamiltonian using Pauli matrices
function generate_Hamiltonian_with_Pauli_matrices(N, J, h)
    Id = [1 0; 0  1]
    Sp = [0 1; 0  0]
    Sm = [0 0; 1  0]
    Sz = [1 0; 0 -1]
    
    # vector of operators: [σᶻ, σᶻ, id, ...]
    fst_term_ops = fill(Id, N)
    fst_term_ops[1] = Sp
    fst_term_ops[2] = Sm
    
    fst_term_conj = fill(Id, N)
    fst_term_conj[1] = Sm
    fst_term_conj[2] = Sp
    
    # vector of operators: [σˣ, id, ...]
    snd_term_ops = fill(Id, N)
    snd_term_ops[1] = Sz
    
    H = zeros(Float64, 2^N, 2^N)
    
    for i in 1:N
        # tensor multiply all operators
        H += J*foldl(kron, fst_term_ops)
        H += J*foldl(kron, fst_term_conj)
        # cyclic shift the operators
        fst_term_ops = circshift(fst_term_ops,1)
        fst_term_conj = circshift(fst_term_conj,1)
    end
    
    for i in 1:N
        H += h*foldl(kron, snd_term_ops)
        snd_term_ops = circshift(snd_term_ops,1)
    end
    
    return -H
end
# Example:
# 8×8 Matrix{Float64}:
#  -3.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0
#  -0.0  -1.0  -1.0  -0.0  -1.0  -0.0  -0.0  -0.0
#  -0.0  -1.0  -1.0  -0.0  -1.0  -0.0  -0.0  -0.0
#  -0.0  -0.0  -0.0   1.0  -0.0  -1.0  -1.0  -0.0
#  -0.0  -1.0  -1.0  -0.0  -1.0  -0.0  -0.0  -0.0
#  -0.0  -0.0  -0.0  -1.0  -0.0   1.0  -1.0  -0.0
#  -0.0  -0.0  -0.0  -1.0  -0.0  -1.0   1.0  -0.0
#  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0   3.0


# Function generating magnetization
function magnetization(state, basis)
    M = 0.
    for (i, bstate) in enumerate(basis)
        bstate_M = 0.
        for spin in bstate
            bstate_M += (state[i]^2 * (isone(spin) ? 1 : -1))/length(bstate)
        end
        M += bstate_M
    end
    return M
end


# Function generating correlation function
function correlation(state, basis, r)
    G = 0.
    for (i, bstate) in enumerate(basis)
        bstate_G = 0.
        for (spin, spin1) in zip(bstate, tuple(circshift(collect(bstate),r)...))
            bstate_G += (state[i]^2 * (isone(spin) ? 1 : -1) * (isone(spin1) ? 1 : -1))/length(bstate)
        end
        G += bstate_G
    end
    return G
end

end
