using SparseArrays
using FastaIO
using DelimitedFiles
using Distances
using Statistics
# using PyPlot
# using PyCall
using Random
using LinearAlgebra
# using Hungarian
using StatsBase
include("./graph_kNN.jl")
include("./GA_SA.jl")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
The function below takes a distance matrix, performs kNN and then anneals according to n_replicates
It's adapted from Carlos Gandarilla-Perez's GA_IPA_DCA_robust - I removed the initial functions to load data
and calculate the distance and removed the final functions that calculated the IPA.
"""

function robust_GA_SA(k::Int64, n_replicates::Int64, M::Int64, N::Int64, FirstSeqSpec::Array{Int64,1}, LastSeqSpec::Array{Int64,1}, MSeqSpec::Array{Int64,1}, IndexSeqSpec::Array{Int64,1}, dij_A::Array{Int64,2}, dij_B::Array{Int64,2}; T_0 =1.0::Float64, alpha = 0.9999::Float64, n_sweep = 40000::Int64)

	#n_replicates is the number of realizations of GA experiment, each time taking different random matchings as starting points.
	#kNN is the number of nearest neighbor in the kNN graph. But, if 0 use Orthology.
	#GA parameters
	#T_0 is the initial temperatures of the exponential schedule.
	#alpha is a factor between 0 and 1 to set the temperature T = T_0 * alpha^n_sweep.
	#n_sweep is the number of sweeps, a sweep is defined as N pairing updates.

	println("Calculating kNN")
	wij_A, wij_B = kNNprop(k, M, dij_A, dij_B)

	#-------------------------------------------------------------------------------------------------------------------------------------------------------
	println("Perform annealing over the number of replicates")
	dataTPprotBepi = Array{Float64, 2}(undef, 4 + M, n_replicates)
	#save the number of sequences in this replicate
	dataTPprotBepi[1, :] = fill(M, n_replicates)
	#compute the mean number of pairs per species
	dataTPprotBepi[2, :] = fill(M/N, n_replicates)
	println("threading")

	Threads.@threads for k in 1:n_replicates
		# println(k)
		dataTPprotBepi[3:4 + M, k] = SA_replicates(T_0, alpha, n_sweep, M, N, FirstSeqSpec, LastSeqSpec, MSeqSpec, IndexSeqSpec, wij_A, wij_B)
	end

	return dataTPprotBepi

end