using Plots
using LsqFit
using StatsBase

function read_creep_file(filename)

    samples = String[]
    stress = Float64[]
    temp = Float64[]
    strain_rate = Float64[]

    open(filename) do io
        for line in eachline(io)

            line = strip(line)

            # ignora righe vuote o header
            if isempty(line) || startswith(line, "Olivine") || startswith(line, "sample")
                continue
            end

            cols = split(line)

            # servono almeno 4 colonne numeriche
            if length(cols) < 4
                continue
            end

            push!(samples, cols[1])
            push!(stress, parse(Float64, cols[2]))
            push!(temp, parse(Float64, cols[3]))
            push!(strain_rate, parse(Float64, cols[4]))
        end
    end

    return samples, stress, temp, strain_rate
end


# --- lettura dati ---
file = "Darot1980_clean.txt"
samples, stress, temp, strain_rate = read_creep_file(file)

R =1 
invT = 1.0 ./ (temp .+ 273.15) # convert in K
log_SR = log10.(strain_rate)

nt = 1000

ns = zeros(length(nt))
Qs = zeros(length(nt))
A0s = zeros(length(nt))


A = zeros(length(stress), 3)
A[:,1] .=1.0
A[:,2] .= log10.(stress)
A[:,3] .= - invT
b = log_SR

x = A\b

n = x[2]
Q = x[3] * R
A0 = 10.0^x[1]





n = mean(ns)
Q = mean(Qs)
A0 = mean(A0s)

σn = std(ns)
σQ = std(Qs)
σA0 = std(A0s)

println("Flow law parameters:")
println("n = $n ± $σn")
println("Q = $Q ± $σQ J/mol")
println("A0 = $A0 ± $σA0 s⁻¹ MPa⁻ⁿ")