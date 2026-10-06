using Test
using ForwardDiff
using ClimFluids

function test()
    params_gas = (kappa = 2 / 7, Cp = 1000, p0 = 1e5, T0 = 300, nu = 0.1)
    params_gas = (kappa0 = params_gas.kappa, Cp0 = params_gas.Cp, params_gas...)
    params_water = (Cp0 = 4000, p0 = 1.1e5, T0 = 300, q0 = 35, v0 = 1e-3, α_p = 1e-9, α_T = 2e-4, α_q = 1e-3, α_TT = 1e-5, γ = 1e-8, μ0 = 0., M_s = 3.14038218e-2, R = 8.31446261815324)

    p, T, q = 1.e5, 295, 34

    # gas:
    IPG = ClimFluids.IdealPerfectGas
    CPV = ClimFluids.VarCpPerfectGas
    
    # liquid:
    NBF = ClimFluids.NonlinearBinaryFluid
    
    consvars_single_gas = Dict( IPG => (:temperature, :entropy, :enthalpy), 
                                CPV => (:temperature,))
    consvars_binary_liquid = Dict( NBF => (:entropy, :potential_temperature))
    
    prec = Float64
    params_gas      = map(prec, params_gas)
    params_water    = map(prec, params_water)
    params_gas      = (; prec, params_gas...)
    params_water    = (; prec, params_water...)

    p, T, q = map(prec, (p, T, q))

    for Fluid in keys(consvars_single_gas)
        for consvar in consvars_single_gas[Fluid]
            ClimFluids.test_fluid(Fluid((; consvar, params_gas...)), p, T)
        end
    end
    for Fluid in keys(consvars_binary_liquid)
        for consvar in consvars_binary_liquid[Fluid]
            ClimFluids.test_fluid(Fluid((; consvar, params_water...)), p, T, q)
        end
    end
end


test()

