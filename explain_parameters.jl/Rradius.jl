using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,BenchmarkTools



αrange=collect(0.5:0.025:1.2)
Rrange=collect(0.5:0.5:14)

@load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[end]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 

p=plot(Rrange,avHY[end,:],grid=false,xticks = ([ 2.5, 5.0, 7.5, 10.0,12.5], [ L"2.5", L"5.0", L"7.5", L"10.0",L"12.5"]),yticks = ([0, 50, 100, 150, 200], [L"0", L"50", L"100", L"150", L"200"]),color=:black,lw=1.5,legend=false,xlabel=L"R",ylabel=L"\mathrm{Average\,\,  degree}", dpi=1200)
savefig(p, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/Rradius.png")
