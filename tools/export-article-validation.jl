# Run from the workspace root with --project=project_link. Julia is needed
# only to refresh validation fixtures, never for installation or R usage.
using LinearAlgebra, Serialization, DelimitedFiles
import RobustSandwichRV
const RS = RobustSandwichRV
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT,"empirical_analysis","delta_omega","bootstrap.jl"))
const D = DeltaOmegaCalibration
const E = PracticalOmegaEnvelope
const O = ObservedOmegaAnalysis
BLAS.set_num_threads(1)
Base.eval(RS, :(struct ProfileOmegaCap
    a0::Float64
    L::Float64
    beta::Float64
    k::Float64
end))
out = joinpath(ROOT, "outputs", "r_package_v05_validation_20261016", "article_inputs")
mkpath(out)
rows = E._read_csv(joinpath(ROOT,"new_article_20261003","generated","table3_casewise_results.csv")).rows
selected = rows
exported = NamedTuple[]
for ds in D.DATASETS
    cases = filter(r -> r["dataset_id"] == ds, selected)
    isempty(cases) && continue
    ctx = D.load_context(ds;root=ROOT)
    columns = ["v" * lpad(i, 3, '0') for i in eachindex(ctx.design)]
    writedlm(joinpath(out, ds * "_input.csv"), vcat(permutedims(vcat(["y"],columns)),hcat(ctx.y,ctx.W)), ',')
    D.write_rows(joinpath(out,ds * "_columns.csv"), [(;name=columns[i],original=ctx.design[i],scalar=ctx.design[i] in ctx.targets) for i in eachindex(columns)])
    for row in cases
        estimator = row["estimator"]; target = row["target_id"]
        rule = ctx.rules[estimator]
        key = ds * "__" * estimator * "__" * replace(target,r"[^A-Za-z0-9_]+"=>"_")
        vcov = startswith(estimator,"CR1") ? "CR1" : estimator
        if vcov == "CR1"
            writedlm(joinpath(out,key * "_cluster.csv"),rule.operator.codes,',')
        elseif vcov == "TWCR1_raw"
            writedlm(joinpath(out,key * "_cluster.csv"),hcat(rule.operator.codes1,rule.operator.codes2),',')
        end
        checkpoint = deserialize(joinpath(ROOT,"new_article_20261003","empirical","profile_calibration",row["checkpoint"]))
        delta = deserialize(joinpath(ROOT,"new_article_20261003","empirical","results_final","checkpoints",
          key * "__delta_omega__scalar.jls"))
        sandrows=filter(r->r["dataset_id"]==ds && r["estimator"]==estimator && r["target_id"]==target && string(r["Kmax"])=="Inf",
          E._read_csv(joinpath(ROOT,"empirical_rvr_geometry_kmax_20261002","runs","full_20261002_211533","results.csv")).rows)
        sand = only(sandrows)
        sandz = deserialize(joinpath(ROOT,"empirical_rvr_geometry_kmax_20261002","runs","full_20261002_211533",sand["witness_file"])).z
        writedlm(joinpath(out,key * "_witness.csv"),hcat(sandz,delta.result.z,checkpoint.result.z),',')
        push!(exported,(;dataset=ds,key,estimator,vcov,target=columns[only(findall(==(target),ctx.design))],
          critical=E._float(row["critical_robust"]),omega=E._float(row["omega_0"]),
          L=E._float(row["log_increase_max"]),beta=E._float(row["slope_max"]),
          profile=E._float(row["profile_lower"]),delta=E._float(row["delta_lower"]),
          delta_upper=E._float(row["delta_upper"]),sand_lower=E._float(row["sand_lower"]),
          sand_upper=E._float(row["sand_upper"])))
    end
    println("EXPORTED ",ds," ",length(cases));flush(stdout)
end
D.write_rows(joinpath(out,"cases.csv"),exported)
println("COMPLETE ",length(exported))
