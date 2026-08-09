using DelimitedFiles
using RobustSandwichRV

package_root = normpath(joinpath(@__DIR__, ".."))
input_path = joinpath(package_root, "inst", "extdata", "julia_reference_input.csv")
output_path = joinpath(package_root, "inst", "extdata", "julia_reference_expected.csv")

raw, header = readdlm(input_path, ',', Float64; header=true)
names = vec(String.(header))
column(name) = raw[:, findfirst(==(name), names)]

y = column("y")
d = column("d")
x1 = column("x1")
x2 = column("x2")
z = column("z")
cluster = Int.(round.(column("cluster")))
X = hcat(ones(length(y)), x1, x2)
cache = prepare_model(y, d, X)

rows = Vector{Vector{Any}}()
for vcov in (:classical, :HC0, :HC1, :CR0, :CR1)
    labels = vcov in (:CR0, :CR1) ? cluster : nothing
    result = evaluate_confounder(cache, z; vcov=vcov, cluster=labels)
    reference = evaluate_confounder_reference(y, d, X, z; vcov=vcov, cluster=labels)
    push!(rows, Any[
        String(vcov), result.tau, result.se_classical, result.se_robust,
        result.omega, result.t_classical, result.t_robust,
        result.r2_d, result.r2_y,
        reference.tau, reference.se_robust, reference.omega,
    ])
end

open(output_path, "w") do io
    println(io, join([
        "vcov", "tau", "se_classical", "se_robust", "omega",
        "t_classical", "t_robust", "r2_d", "r2_y",
        "reference_tau", "reference_se_robust", "reference_omega",
    ], ','))
    for row in rows
        println(io, join(map(value -> value isa AbstractString ? value :
            string(Float64(value)), row), ','))
    end
end
