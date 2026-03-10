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

model(x,p) = p[1] .+ p[2] .* x   # p[1] = intercetta, p[2] = slope


# --- lettura dati ---
file = "Darot1980_clean.txt"
samples, stress, temp, strain_rate = read_creep_file(file)

R =1 
invT = 1.0 ./ (temp .+ 273.15) # convert in K
log_SR = log10.(strain_rate)


ns = zeros(length(stress))
Qs = zeros(length(stress))
A1s = zeros(length(stress))
A2s = zeros(length(stress))
A0s = zeros(length(stress))

unique_samples = unique(samples)
unique_temp = unique(temp)
sample_idx = Dict(s => findall(==(s), samples) for s in unique_samples)
sample_idx_T = Dict(t => findall(==(t), temp) for t in unique_temp)

# marker diversi per ogni sample
marker_list = [:circle, :square, :diamond, :utriangle, :star5]
sample_marker = Dict(unique_samples[i] => marker_list[i] for i in eachindex(unique_samples))

# colori diversi per temperatura
palette_colors = palette(:viridis, length(unique_temp))
temp_color = Dict(unique_temp[i] => palette_colors[i] for i in eachindex(unique_temp))


# --- plot ---
p = plot(
    xlabel = "Strain rate (s⁻¹)",
    ylabel = "Deviatoric stress (MPa)",
    xscale = :log10,
    legend = :outerright
)

for i in eachindex(stress)

    scatter!(
        p,
        [strain_rate[i]],
        [stress[i]],
        marker = sample_marker[samples[i]],
        color = temp_color[temp[i]],
        label = "$(samples[i]), $(temp[i])°C"
    )

end

#display(p)


# log ̇ε = log A2 + Q/RT
p1 = plot(layout=(2,2), size=(1600,1200))


for (i,s) in enumerate(unique_samples)

    idx = sample_idx[s]

    fit = curve_fit(model, invT[idx], log_strain[idx], [0.0,1.0])

    intercept = fit.param[1]
    A2s[idx] .= 10.0^intercept
    slope = fit.param[2]
    Qs[idx] .= -slope * R


    # plot 

    scatter!(
        p1[i],
        invT[idx],
        log_strain[idx],
        yerror = res[idx],  
        xlabel = "1/T (1/K)",
        ylabel = "log₁₀ strain rate (s⁻¹)",
        title = "$s",
        legend = false
    )


    xs = range(minimum(invT[idx]), maximum(invT[idx]), length=10)
    ys = intercept .+ slope .* xs


    plot!(p1[i], xs, ys)

end

#display(p1)

res = zeros(length(stress))

# log ̇ε = log A1 + n log σ

p2 = plot(layout=(2,3), size=(1600,1200))

for (i,s) in enumerate(unique_temp)

    idx = sample_idx_T[s]


    fit = curve_fit(model, log10.(stress[idx]), log_strain[idx], [0.0,1.0])

    intercept = fit.param[1]
    A1s[idx] .= 10.0^intercept
    slope = fit.param[2]
    ns[idx] .= slope


    #plot
    scatter!(
        p2[i],
        log10.(stress[idx]),
        log_strain[idx],
        yerror = res[idx],
        xlabel = "log₁₀ Deviatoric stress (MPa)",
        ylabel = "log₁₀ strain rate (s⁻¹)",
        title = "$s",
        legend = false
        )

    xs = range(minimum(log10.(stress[idx])), maximum(log10.(stress[idx])), length=10)
    ys = intercept .+ slope .* xs
    plot!(p2[i], xs, ys)

end

#display(p2)

n = mean(ns)
Q = mean(Qs)
A1 = mean(A1s)
A2 = mean(A2s)

A0s .= strain_rate ./ (stress.^ns .* exp.(-Qs ./ (R .* (temp .+ 273.15))))
A0 = mean(A0s)

σn = std(ns)
σQ = std(Qs)
σA0 = std(A0s)

println("Flow law parameters:")
println("n = $n ± $σn")
println("Q = $Q ± $σQ J/mol")
println("A0 = $A0 ± $σA0 s⁻¹ MPa⁻ⁿ")