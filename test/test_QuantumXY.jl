using Test
using QuantumXY

@testset "QuantumXY.jl" begin
	@testset "Correct basis" begin
		# Function generating full basis
		@test generate_full_basis(3) == [(1, 1, 1),(0, 1, 1),(1, 0, 1),(0, 0, 1),(1, 1, 0),(0, 1, 0),(1, 0, 0),(0, 0, 0)]
	end
	@testset "Equivalence of Hamiltonians" begin
		# Equivalence of Hamiltonians
		@test isapprox(generate_Hamiltonian(generate_full_basis(3), 1, 1), generate_Hamiltonian_with_Pauli_matrices(3,1,1))
	end
	@testset "Correct magnetization" begin
		# Magnetization
		@test magnetization([1,0,0,0,0,0,0,0], generate_full_basis(3)) == 1.0
	end
	@testset "Corrected correlation" begin
		# Connected correlation
		@test correlation([1,0,0,0,0,0,0,0], generate_full_basis(3), 1) == 1.0
	end
end
