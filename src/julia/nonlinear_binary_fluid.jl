"""
    fluid = NonlinearBinaryFluid(consvar, params)
    fluid = NonlinearBinaryFluid(options)
Return an object `fluid` describing the thermodynamics of a nonlinear binary fluid.
The fluid includes cabbeling and thermobaric terms.
`fluid` can then be used as first argument of thermodynamic functions.
Argument `consvar` defines the conservative variable chosen. Valid values are `:entropy`, `:potential_temperature`.

Named tuple `params` contains the reference values `p0` (pressure), `T0` (temperature), `q0` (composition, kg/kg),
`v0` (specific volume) and `Cp0` (specific heat capacity at `p0`); the thermal expansion `α_T`, haline contraction `α_q`,
compressibility `α_p`, cabbeling `α_TT` and thermobaric `γ` coefficients; the gas constant `R` and the mean molar
mass of the solute `M_s`. Specific volume is

    v = v0 (1 + α_T (1 + γ (p - p0)) (θ - T0) - α_q (q - q0) - α_p (p - p0) + α_TT (θ - T0)² / 2)

and specific entropy s = Cp0 log(θ / T0) + s_mix(q), with s_mix(q) = -(R / M_s) q (log q - 1).

Arguments can also be passed as a named tuple `options`, including `consvar`.
"""
struct NonlinearBinaryFluid{CV, F} <: BinaryFluid{F}
    p0::F
    T0::F
    q0::F
    v0::F
    Cp0::F
    α_T::F
    α_q::F
    α_p::F
    α_TT::F
    γ::F
    M_s::F
    R::F
    R_s::F
    function NonlinearBinaryFluid(consvar::Symbol, (; p0, T0, q0, v0, Cp0, α_T, α_q, α_p, α_TT, γ, R, M_s))
        @assert consvar in (:entropy, :potential_temperature)
        new{consvar, typeof(p0)}(p0, T0, q0, v0, Cp0, α_T, α_q, α_p, α_TT, γ, M_s, R, R/M_s)
    end
end
NonlinearBinaryFluid(params) = NonlinearBinaryFluid(params.consvar, params)

const NBF   = NonlinearBinaryFluid
const NBFS  = NonlinearBinaryFluid{:entropy}
const NBFP  = NonlinearBinaryFluid{:potential_temperature}

const PVCons = NamedTuple{(:p, :v, :consvar)}

@fastmath @muladd @inlineall begin
    ## private functions
    xlnx(x, prec) = x * log(x + eps(prec))

    # p(v, T, q) from solving v(p, θ(p, T), q) - v = 0
    function NBF_p_vTq(fluid, v, T, q)
        (; p0) = fluid
        p = p0
        for _ in 1:4
            theta = NBF_θ_pT(fluid, p, T)
            f = NBF_v_pθq(fluid, p, theta, q) - v
            fprime = NBF_∂v_∂p(fluid, theta) + NBF_∂v_∂θ(fluid, p, theta) * NBF_∂θ_∂p(fluid, p, T, theta)
            p = p - f / fprime
        end
        return p
    end
    # p(v, θ, q)
    NBF_p_vθq((; p0, α_T, γ, T0, α_p, v0, α_q, q0, α_TT), v, θ, q) = p0 + inv( α_T * γ * (θ - T0) - α_p ) * ( (v - v0)/v0 + α_q * (q - q0) - α_T * (θ - T0) - 0.5 * α_TT * (θ - T0)^2 )
    # T(p, θ)
    NBF_T_pθ((; v0, Cp0, α_T, γ, α_TT, p0, T0), p, θ) = (θ / Cp0) * NBF_exner_pθ((; Cp0, v0, α_T, γ, p0, α_TT, T0), p, θ)
    # ∂T/∂θ(p, θ)
    NBF_∂T_∂θ((; v0, Cp0, α_T, γ, α_TT, p0, T0), p, θ) = 1 + v0 * ( α_T * (1 + 0.5 * γ * (p - p0)) + α_TT * (2*θ - T0) ) * (p - p0) / Cp0
    # θ(p, T): T = θ (1 + a + b (θ - T0)) is quadratic in θ, b θ² + B θ - T = 0 with B = 1 + a - b T0,
    # whose root θ = 2T / (B + √(B² + 4bT)) reduces to θ = T / B as b → 0
    function NBF_θ_pT((; v0, Cp0, α_T, γ, α_TT, p0, T0), p, T)
        a = v0 * α_T * (1 + 0.5 * γ * (p - p0)) * (p - p0) / Cp0
        b = v0 * α_TT * (p - p0) / Cp0
        B = 1 + a - b * T0
        return 2T / (B + sqrt(B^2 + 4b * T))
    end
    # ∂θ/∂p(p, T) at fixed T, = -(∂T/∂p)_θ / (∂T/∂θ)_p
    NBF_∂θ_∂p((; v0, Cp0, α_T, γ, α_TT, p0, T0), p, T, theta) = - ( (v0/Cp0) * ( α_T + α_T * γ * (p - p0) + α_TT * (theta - T0) ) * theta) / ( 1 + (v0 / Cp0) * ( α_T + 0.5 * α_T * γ * (p - p0) + α_TT * (2*theta - T0) ) * (p - p0) )
    # θ(p, v, q): with x = θ - T0, v(p, θ, q) = v gives ½ α_TT x² + A x - R = 0, A = α_T (1 + γ (p - p0)),
    # R = v/v0 - 1 + α_q (q - q0) + α_p (p - p0); root x = 2R / (A + √(A² + 2 α_TT R)) (warm side of any density maximum)
    function NBF_θ_pvq((; p0, T0, q0, v0, α_p, α_T, α_q, γ, α_TT), p, v, q)
        A = α_T * (1 + γ * (p - p0))
        R = v / v0 - 1 + α_q * (q - q0) + α_p * (p - p0)
        return T0 + 2R / (A + sqrt(A^2 + 2α_TT * R))
    end
    # θ(s, q)
    NBF_θ_sq(fluid, s, q) = fluid.T0 * exp( (s - NBF_s_mix(fluid, q)) / fluid.Cp0)
    # ∂θ/∂s(s, q)
    NBF_∂θ_∂s(fluid, s, q) = NBF_θ_sq(fluid, s, q) / fluid.Cp0
    # ∂θ/∂q(s, q)
    NBF_∂θ_∂q(fluid, s, q) = - NBF_∂s_mix_∂q(fluid, q) * NBF_θ_sq(fluid, s, q) / fluid.Cp0
    # s(θ, q)
    NBF_s_θq(fluid, θ, q) = NBF_s_mix(fluid, q) + fluid.Cp0 * log(θ / fluid.T0)
    # s_mix(q) (q in kg/kg)
    function NBF_s_mix(fluid::NonlinearBinaryFluid{CV, F}, q) where {CV, F}
        (; R_s) = fluid
        return - R_s * ( xlnx(q, F) - q )
    end
    # s_mix'(q)
    NBF_∂s_mix_∂q((; R_s), q) = - R_s * log(q )
    # s_mix''(q)
    NBF_∂s_mix_∂qq((; R_s), q) = - R_s / q
    # v(p, θ, q) = ∂h/∂p(p, θ, q)
    NBF_v_pθq((; p0, T0, q0, v0, α_p, α_T, α_q, γ, α_TT), p, θ, q) = v0 * ( 1 + α_T * (1 + γ * (p - p0)) * (θ - T0) - α_q * (q - q0) - α_p * (p - p0) + 0.5 * α_TT * (θ - T0)^2 )
    # ∂v/∂p(p, θ, q)
    NBF_∂v_∂p((; v0, α_p, α_T, γ, T0), θ) = - v0 * α_p + v0 * α_T * γ * (θ - T0)
    # ∂v/∂θ(p, θ, q)
    NBF_∂v_∂θ((; v0, α_T, γ, α_TT, p0, T0), p, θ) = v0 * α_T * ( 1 + γ * (p - p0)) + v0 * α_TT * (θ - T0)
    # π = Cp0.T / θ = ∂h/∂θ(p, θ, q)
    NBF_exner_pθ((; Cp0, v0, α_T, γ, p0, α_TT, T0), p, θ) = Cp0 + v0 * ( α_T * (1 + 0.5 * γ * (p - p0)) + α_TT * (θ - T0) ) * (p - p0)
    # h(p, θ, q)
    NBF_h_pθq((; p0, T0, q0, v0, α_p, α_T, α_q, γ, α_TT, Cp0), p, θ, q) = Cp0 * (θ - T0) + v0 * ( 1 + α_T * (1 + 0.5 * γ * (p - p0)) * (θ - T0) - α_q * (q - q0) - 0.5 * α_p * (p - p0) + 0.5 * α_TT * (θ - T0)^2 ) * (p - p0)
    # c2(p, θ, q)
    NBF_c2_pθq((; p0, T0, q0, v0, α_p, α_T, α_q, γ, α_TT), p, θ, q) = - NBF_v_pθq((; p0, T0, q0, v0, α_p, α_T, α_q, γ, α_TT), p, θ, q)^2 / NBF_∂v_∂p((; v0, α_p, α_T, γ, T0), θ)
    # μ(p, θ, q) = ∂h/∂q(p, s, q) = ∂h/∂q(p, θ, q) + ∂h/∂θ(p, θ, q) ∂θ/∂q(s, q)
    NBF_chemical_potential_pθq((; v0, p0, T0, α_T, α_q, γ, α_TT, R_s, Cp0), p, θ, q) = - v0 * α_q * (p - p0) - NBF_∂s_mix_∂q((; R_s), q) * NBF_T_pθ((; v0, Cp0, α_T, γ, α_TT, p0, T0), p, θ)
    # μ(p, T, q) = - v0 α_q (p - p0) - s_mix'(q) T, which needs no inversion for θ
    NBF_chemical_potential_pTq((; v0, p0, α_q, R_s), p, T, q) = - v0 * α_q * (p - p0) - NBF_∂s_mix_∂q((; R_s), q) * T
    # Cp(p, θ) = T (∂s/∂T)_p = π / (∂T/∂θ)_p, since T = θ π / Cp0 and s_θ = Cp0 / θ; equals Cp0 at p = p0
    NBF_Cp_pθ(fluid, p, θ) = NBF_exner_pθ(fluid, p, θ) / NBF_∂T_∂θ(fluid, p, θ)
    # μθ(p, θ, q) = ∂h/∂q(p, θ, q)
    NBF_modified_chemical_potential_pq((; v0, α_q, p0), p, q) = - v0 * α_q * (p - p0) * one(q)
    # q(p, v, θ) 
    NBF_composition_pvθ((; p0, T0, q0, v0, α_p, α_T, α_q, γ, α_TT), p, v, θ) = q0 + (-(v - v0)/v0 + α_T * (1 + γ * (p - p0)) * (θ - T0) - α_p * (p - p0) + 0.5 * α_TT * (θ - T0)^2 )/α_q
    # q(p, v, s) from solving v(p, θ(s, q), q) - v = 0, since θ(s, q) depends on q through s_mix(q).
    # v - v(p, θ(s, q), q) is monotone in q, so the root is unique; steps are limited to q/2 to keep q > 0.
    function NBF_composition_pvs(fluid::NonlinearBinaryFluid{CV, F}, p, v, s) where {CV, F}
        q = fluid.q0
        for _ in 1:8
            θ = NBF_θ_sq(fluid, s, q)
            f = NBF_v_pθq(fluid, p, θ, q) - v
            fprime = NBF_∂v_∂θ(fluid, p, θ) * NBF_∂θ_∂q(fluid, s, q) - fluid.v0 * fluid.α_q
            q = max(q - f / fprime, q / 2)
        end
        return q
    end

    ## all consvar
    # defined throughout ClimFluids:
    specific_entropy(           fluid::NBF, (p, T, q)::PTQ)     = NBF_s_θq(fluid, NBF_θ_pT(fluid, p, T), q)
    specific_enthalpy(          fluid::NBF, (p, T, q)::PTQ)     = NBF_h_pθq(fluid, p, NBF_θ_pT(fluid, p, T), q)
    specific_internal_energy(   fluid::NBF, (p, T, q)::PTQ)     = NBF_h_pθq(fluid, p, NBF_θ_pT(fluid, p, T), q) - p * NBF_v_pθq(fluid, p, NBF_θ_pT(fluid, p, T), q)
    potential_temperature(      fluid::NBF, (p, T, q)::PTQ)     = NBF_θ_pT(fluid, p, T)
    potential_enthalpy(         fluid::NBF, (p, T, q)::PTQ)     = NBF_h_pθq(fluid, fluid.p0, NBF_θ_pT(fluid, p, T), q)
    potential_volume(           fluid::NBF, (p, T, q)::PTQ)     = NBF_v_pθq(fluid, fluid.p0, NBF_θ_pT(fluid, p, T), q)
    heat_capacity(              fluid::NBF, (p, T, q)::PTQ)     = NBF_Cp_pθ(fluid, p, NBF_θ_pT(fluid, p, T))
    pressure(                   fluid::NBF, (v, T, q)::VTQ)     = NBF_p_vTq(fluid, v, T, q)
    sound_speed2(               fluid::NBF, (p, T, q)::PTQ)     = NBF_c2_pθq(fluid, p, NBF_θ_pT(fluid, p, T), q)
    
    # defined here only:
    potential_temperature(      fluid::NBF, (v, T, q)::VTQ)     = NBF_θ_pvq(fluid, NBF_p_vTq(fluid, v, T, q), v, q)
    specific_volume(            fluid::NBF, (p, T, q)::PTQ)     = NBF_v_pθq(fluid, p, NBF_θ_pT(fluid, p, T), q)
    chemical_potential(         fluid::NBF, (p, T, q)::PTQ)     = NBF_chemical_potential_pTq(fluid, p, T, q)
    specific_entropy(           fluid::NBF, (p, θ, q)::PThQ)    = NBF_s_θq(fluid, θ, q)

    function chemical_potential_derivatives(fluid::NBF, (p, T, q)::PTQ) # returns μ(p, T, q), ∂μ/∂p(p, T, q), ∂μ/∂T(p, T, q), ∂μ/∂q(p, T, q)
        μ = NBF_chemical_potential_pTq(fluid, p, T, q)
        ∂μ_∂p = - fluid.v0 * fluid.α_q
        ∂μ_∂T = - NBF_∂s_mix_∂q(fluid, q)
        ∂μ_∂q = - NBF_∂s_mix_∂qq(fluid, q) * T
        return μ, ∂μ_∂p, ∂μ_∂T, ∂μ_∂q
    end

    ## consvar = :entropy
    # PTQ
    conservative_variable(      fluid::NBFS, (p, T, q)::PTQ)    = NBF_s_θq(fluid, NBF_θ_pT(fluid, p, T), q)
    conjugate_variable(         fluid::NBFS, (p, T, q)::PTQ)    = T
    dcons_dq(                   fluid::NBFS, (p, T, q)::PTQ)    = NBF_∂s_mix_∂q(fluid, q) # with p, T const
    modified_chemical_potential(fluid::NBFS, (p, T, q)::PTQ)    = NBF_chemical_potential_pTq(fluid, p, T, q)

    # PConsQ
    temperature(                fluid::NBFS, (p, s, q)::PConsQ) = NBF_T_pθ(fluid, p, NBF_θ_sq(fluid, s, q))
    specific_enthalpy(          fluid::NBFS, (p, s, q)::PConsQ) = NBF_h_pθq(fluid, p, NBF_θ_sq(fluid, s, q), q)
    specific_volume(            fluid::NBFS, (p, s, q)::PConsQ) = NBF_v_pθq(fluid, p, NBF_θ_sq(fluid, s, q), q)
    heat_capacity(              fluid::NBFS, (p, s, q)::PConsQ) = NBF_Cp_pθ(fluid, p, NBF_θ_sq(fluid, s, q))
    modified_chemical_potential(fluid::NBFS, (p, s, q)::PConsQ) = NBF_chemical_potential_pθq(fluid, p, NBF_θ_sq(fluid, s, q), q)
    dcons_dq(                   fluid::NBFS, (p, s, q)::PConsQ) = NBF_∂s_mix_∂q(fluid, q) # with p, T const
    potential_temperature(      fluid::NBFS, (p, s, q)::PConsQ) = NBF_θ_sq(fluid, s, q)
    conjugate_variable(         fluid::NBFS, (p, s, q)::PConsQ) = NBF_T_pθ(fluid, p, NBF_θ_sq(fluid, s, q))

    function exner_functions(fluid::NBFS, (p, s, q)::PConsQ) # returns h, v, conjvar = T
        θ = NBF_θ_sq(fluid, s, q)
        h = NBF_h_pθq(fluid, p, θ, q)
        v = NBF_v_pθq(fluid, p, θ, q)
        conjvar = NBF_T_pθ(fluid, p, θ)
        return h, v, conjvar
    end
    function volume_functions(fluid::NBFS, (p, s, q)::PConsQ) # returns v, ∂v/∂p, ∂v/∂s, ∂v/∂q i.e. h_pp, h_ps, h_pq
        θ = NBF_θ_sq(fluid, s, q)
        v = NBF_v_pθq(fluid, p, θ, q)
        dv_dp = NBF_∂v_∂p(fluid, θ)
        dv_ds = (θ / fluid.Cp0) * NBF_∂v_∂θ(fluid, p, θ)
        dv_dq = -fluid.v0 * fluid.α_q + NBF_∂v_∂θ(fluid, p, θ) * NBF_∂θ_∂q(fluid, s, q)
        return v, dv_dp, dv_ds, dv_dq
    end

    # VCons
    temperature(    fluid::NBFS, (v, s, q)::VConsQ)             = NBF_T_pθ(fluid, NBF_p_vθq(fluid, v, NBF_θ_sq(fluid, s, q), q), NBF_θ_sq(fluid, s, q))
    pressure(       fluid::NBFS, (v, s, q)::VConsQ)             = NBF_p_vθq(fluid, v, NBF_θ_sq(fluid, s, q), q)

    ## consvar = :potential_temperature
    # PTQ
    conservative_variable(  fluid::NBFP, (p, T, q)::PTQ)        = NBF_θ_pT(fluid, p, T)
    conjugate_variable(     fluid::NBFP, (p, T, q)::PTQ)        = NBF_exner_pθ(fluid, p, NBF_θ_pT(fluid, p, T))
    modified_chemical_potential(fluid::NBFP, (p, T, q)::PTQ)    = NBF_modified_chemical_potential_pq(fluid, p, q)
    dcons_dq(               fluid::NBFP, (p, T, q)::PTQ)        = zero(T) # with p, T const

    # PConsQ
    temperature(            fluid::NBFP, (p, θ, q)::PConsQ)         = NBF_T_pθ(fluid, p, θ)
    specific_enthalpy(      fluid::NBFP, (p, θ, q)::PConsQ)         = NBF_h_pθq(fluid, p, θ, q)
    specific_volume(        fluid::NBFP, (p, θ, q)::PConsQ)         = NBF_v_pθq(fluid, p, θ, q)
    heat_capacity(          fluid::NBFP, (p, θ, q)::PConsQ)         = NBF_Cp_pθ(fluid, p, θ)
    modified_chemical_potential(fluid::NBFP, (p, θ, q)::PConsQ)     = NBF_modified_chemical_potential_pq(fluid, p, q)
    dcons_dq(               fluid::NBFP, (p, θ, q)::PConsQ)         = zero(θ) # with p, T const
    potential_temperature(  fluid::NBFP, (p, θ, q)::PConsQ)         = θ
    conjugate_variable(     fluid::NBFP, (p, θ, q)::PConsQ)         = NBF_exner_pθ(fluid, p, θ)

    function exner_functions(fluid::NBFP, (p, θ, q)::PConsQ) # returns h, v, conjvar = Cp0 * T / θ
        h = NBF_h_pθq(fluid, p, θ, q)
        v = NBF_v_pθq(fluid, p, θ, q)
        conjvar = NBF_exner_pθ(fluid, p, θ)
        return h, v, conjvar
    end
    function volume_functions(fluid::NBFP, (p, θ, q)::PConsQ) # returns v(p, θ, q), ∂v/∂p, ∂v/∂θ, ∂v/∂q
        v = NBF_v_pθq(fluid, p, θ, q)
        ∂v_∂p = NBF_∂v_∂p(fluid, θ)
        ∂v_∂θ = NBF_∂v_∂θ(fluid, p, θ)
        ∂v_∂q = - fluid.v0 * fluid.α_q
        return v, ∂v_∂p, ∂v_∂θ, ∂v_∂q
    end

    # VCons
    temperature(    fluid::NBFP, (v, θ, q)::VConsQ)             = NBF_T_pθ(fluid, NBF_p_vθq(fluid, v, θ, q), θ)
    pressure(       fluid::NBFP, (v, θ, q)::VConsQ)             = NBF_p_vθq(fluid, v, θ, q)

    ## Fallback: convert to (p, T, q) if not implemented above
    canonical_state(fluid::NBF, (p, T, q)::PTQ)                 = (p = p, T = T, q = q)
    canonical_state(fluid::NBF, (p, consvar, q)::PConsQ)        = (p, T = NBF_T_pθ(fluid, p, potential_temperature(fluid, (; p, consvar, q))), q)
    canonical_state(fluid::NBF, (p, s, q)::PSQ)                 = (p, T = NBF_T_pθ(fluid, p, NBF_θ_sq(fluid, s, q)), q)
    canonical_state(fluid::NBF, (p, θ, q)::PThQ)                = (p, T = NBF_T_pθ(fluid, p, θ), q)
    canonical_state(fluid::NBF, (p, v, q)::PVQ)                 = (p, T = NBF_T_pθ(fluid, p, NBF_θ_pvq(fluid, p, v, q)), q)
    canonical_state(fluid::NBF, (v, T, q)::VTQ)                 = (p = pressure(fluid, (; v, T, q)), T, q)
    canonical_state(fluid::NBF, (v, s, q)::VSQ)                 = canonical_state_NBF_vθ(fluid, (v, NBF_θ_sq(fluid, s, q), q))
    canonical_state(fluid::NBFS,(v, s, q)::VConsQ)              = canonical_state_NBF_vθ(fluid, (v, NBF_θ_sq(fluid, s, q), q))
    canonical_state(fluid::NBFP,(v, θ, q)::VConsQ)              = canonical_state_NBF_vθ(fluid, (v, θ, q))

    function canonical_state_NBF_vθ(fluid, (v, θ, q))
        p = NBF_p_vθq(fluid, v, θ, q)
        T = NBF_T_pθ(fluid, p, θ)
        return (; p, T, q)
    end

    exner_functions(fluid::NBF, (p, T, q)::PTQ)     = exner_functions(fluid,   (; p, consvar = conservative_variable(fluid, (; p, T, q)), q))
    volume_functions(fluid::NBF, (p, T, q)::PTQ)    = volume_functions(fluid,  (; p, consvar = conservative_variable(fluid, (; p, T, q)), q))

    # check if needed:
    ds_dq(      fluid::NBF, (p, T, q)::PTQ)     = NBF_∂s_mix_∂q(fluid, q)   # ∂s/∂q with p, T const.
    composition(fluid::NBFS, (p, v, s)::PVCons) = NBF_composition_pvs(fluid, p, v, s)
    composition(fluid::NBFP, (p, v, θ)::PVCons) = NBF_composition_pvθ(fluid, p, v, θ)
end
