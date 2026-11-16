using GeneralizedGauss
using LinearAlgebra

const BENCH_BITS = parse(Int, get(ENV, "GG_BENCH_BITS", "512"))
const BENCH_TOLERANCE = get(ENV, "GG_BENCH_TOLERANCE", "1e-40")
const BENCH_SCALE = get(ENV, "GG_BENCH_SCALE", "none")
const BENCH_ENDPOINT = get(ENV, "GG_BENCH_ENDPOINT", "natural")

# Prescribed moments: no call to compute_moments and no assumed target density.
function fixture(name)
    orders = if name == "mixed23"
        sort(unique([2i + 3j for i in 0:5 for j in 0:5]))
    elseif name == "high201"
        [0, 100, 101, 201]
    elseif name in ("mixed50", "theory50")
        [0, 10, 40, 50]
    elseif name == "radau"
        [0, 2, 5]
    else
        error("Unknown fixture")
    end
    if name == "theory50"
        a = BigFloat(891)/100
        m = [prod((a+i)/(BigFloat(9)+i) for i in 0:k-1; init=BigFloat(1)) for k in orders]
        return orders, m, BigFloat[], BigFloat[]
    elseif name == "radau"
        m = [inv(BigFloat(k + 1)) for k in orders]
        z = (m[3] / m[2])^(inv(BigFloat(3)))
        w = m[2] / z^2
        return orders, m, BigFloat[0,z], BigFloat[1-w,w]
    end
    n = length(orders) ÷ 2
    x = BigFloat.(1:n) ./ (n + 1)
    w = BigFloat.(1:n) ./ sum(1:n)
    m = [sum(w .* x.^k) for k in orders]
    orders, m, x, w
end

function run_case(name, orth, strict)
    orders, m, xref, wref = fixture(name)
    scale = BENCH_SCALE == "first" ? orders[2] : 1
    powers = scale == 1 ? orders : BigFloat.(orders) ./ scale
    funs = Function[let k=k; x -> x^k end for k in powers]
    derivs = Function[let k=k; x -> k == 0 ? zero(x) : k*x^(k-1) end for k in powers]
    basis = quadbasis(funs, derivs, BigFloat(0), BigFloat(1))
    started = time_ns()
    usedbasis, usedmoments = if orth
        qb, transform = orthogonalize_basis(basis;
            quad_order=max(2length(orders), maximum(orders)+1))
        qb, transform*m
    else
        basis, m
    end
    setup_seconds = (time_ns()-started)/1e9
    for side in (:lower, :upper)
        started = time_ns()
        println("START\t",name,"\t",orth,"\t",strict,"\t",side); flush(stdout)
        try
            w,x = compute_gauss_rule(usedbasis, usedmoments;
                principal=side, solver_tolerance=BigFloat(BENCH_TOLERANCE),
                add_endpoint=iseven(length(orders)) && BENCH_ENDPOINT != "natural" ?
                    Symbol(BENCH_ENDPOINT) : nothing,
                intermediate_tolerance=strict ? :strict : nothing,
                lost_digits=0)
            seconds = (time_ns()-started)/1e9
            x = x .^ inv(BigFloat(scale))
            ix = sortperm(x); x=x[ix]; w=w[ix]
            residuals = [abs(sum(w .* x.^k)-m[j]) for (j,k) in enumerate(orders)]
            absres = maximum(residuals)
            relres = maximum(residuals ./ max.(abs.(m), eps(BigFloat)))
            referr = side == :lower && length(x)==length(xref) ?
                max(maximum(abs.(x-xref)), maximum(abs.(w-wref))) : BigFloat(NaN)
            valid = all(isfinite,w) && all(isfinite,x) && all(w .> 0) &&
                all((0 .<= x) .& (x .<= 1)) && all(diff(x) .> 0)
            if iseven(length(orders))
                valid &= side == :lower ?
                    length(x)==length(orders)÷2 && first(x)>0 && last(x)<1 :
                    length(x)==length(orders)÷2+1 && first(x)==0 && last(x)==1
            else
                valid &= length(x)==(length(orders)+1)÷2 &&
                    (side == :lower ? first(x)==0 && last(x)<1 :
                        first(x)>0 && last(x)==1)
            end
            screened = valid && absres <= BigFloat("1e-9") &&
                relres <= BigFloat("1e-10") && (isnan(referr) || referr <= BigFloat("1e-8"))
            println("RESULT\t",name,"\t",orth,"\t",strict,"\t",side,"\t",seconds,
                "\tsetup=",setup_seconds,"\tabsres=",Float64(absres),"\trelres=",Float64(relres),
                "\treferr=",Float64(referr),"\tvalid=",valid,"\tscreened=",screened,
                "\tprevalence=",Float64(1-sum(w.*x)),"\tnodes=",length(x))
            println("RULE\t",name,"\t",orth,"\t",strict,"\t",side,
                "\t",join(string.(x),','),"\t",join(string.(w),','))
        catch e
            println("FAIL\t",name,"\t",orth,"\t",strict,"\t",side,"\t",
                (time_ns()-started)/1e9,"\t",sprint(showerror,e))
        end
        flush(stdout)
    end
end

setprecision(BigFloat,BENCH_BITS) do
    println("CONFIG\tbits=",BENCH_BITS,"\ttolerance=",BENCH_TOLERANCE,
        "\tscale=",BENCH_SCALE,"\tendpoint=",BENCH_ENDPOINT)
    # Compile the main analytic solver paths before reporting timings.
    run_case("radau", false, true)
    for config in ARGS
        name, orth, strict = split(config, ':')
        run_case(name, orth=="qr", strict=="strict")
    end
end
