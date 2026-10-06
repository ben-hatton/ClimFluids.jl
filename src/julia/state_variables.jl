const PV = NamedTuple{(:p, :v)}
const PS = NamedTuple{(:p, :s)}
const PT = NamedTuple{(:p, :T)}
const PTh = NamedTuple{(:p, :theta)}
const PCons = NamedTuple{(:p, :consvar)}
const VT = NamedTuple{(:v, :T)}
const VS = NamedTuple{(:v, :s)}
const VCons = NamedTuple{(:v, :consvar)}

const PVQ = NamedTuple{(:p, :v, :q)}
const PSQ = NamedTuple{(:p, :s, :q)}
const PTQ = NamedTuple{(:p, :T, :q)}
const PThQ = NamedTuple{(:p, :theta, :q)}
const PConsQ = NamedTuple{(:p, :consvar, :q)}
const VTQ = NamedTuple{(:v, :T, :q)}
const VSQ = NamedTuple{(:v, :s, :q)}
const VConsQ = NamedTuple{(:v, :consvar, :q)}

@inlineall begin

    state_variables(::AbstractFluid) = (
        (:p, :v, :q),
        (:p, :T, :q),
        (:p, :s, :q),
        (:p, :consvar, :q),
        (:p, :theta, :q),
        (:v, :T, :q),
        (:v, :s, :q),
        (:v, :consvar, :q),
    )

    # simple fluids do not have composition `q`.
    state_variables(::SimpleFluid) = (
        (:p, :v),
        (:p, :T),
        (:p, :s), # keep ?
        (:p, :consvar),
        (:v, :T),
        (:v, :s),
        (:v, :consvar),
    )

    all_state_functions() = (;
        temperature,
        pressure,
        conservative_variable,
        conjugate_variable,
        specific_volume,
        specific_entropy,
        specific_enthalpy,
        specific_internal_energy,
        potential_volume,
        potential_enthalpy,
        potential_temperature,
        exner_functions,
        volume_functions,
        sound_speed2,
        sound_speed
    )

    binary_state_functions() = (;
        temperature,
        pressure,
        conservative_variable,
        conjugate_variable,
        specific_volume,
        specific_entropy,
        specific_enthalpy,
        specific_internal_energy,
        potential_volume,
        potential_enthalpy,
        potential_temperature,
        exner_functions,
        volume_functions,
        sound_speed2,
        sound_speed,
        heat_capacity,
        modified_chemical_potential,
        dcons_dq,
        chemical_potential_derivatives
    )

    # one could imagine to specialize this function to restrict the list
    state_functions(::AbstractFluid) = all_state_functions()
    state_functions(::BinaryFluid) = binary_state_functions()

    function fallback(fun, fluid::AbstractFluid, state::NamedTuple)
        cstate = canonical_state(fluid, state)
        if cstate == state
            error("Fluid $(typeof(fluid)) does not implement $fun !")
        else
            return fun(fluid, cstate)
        end
    end

    # trivial thermodynamic functions
    temperature(::SimpleFluid, (p,T)::PT) = T
    temperature(::SimpleFluid, (v,T)::VT) = T
    pressure(::SimpleFluid, (p,T)::PT) = p

    temperature(::BinaryFluid, (p,T,q)::PTQ) = T
    temperature(::BinaryFluid, (v,T,q)::VTQ) = T
    pressure(::BinaryFluid, (p,T,q)::PTQ) = p
end

# fallback, only sound_speed2 needs to be implemented
@inline sound_speed(fluid::AbstractFluid, state::NamedTuple) = @fastmath sqrt(sound_speed2(fluid, state))

for fun in propertynames(all_state_functions())
    if fun != :sound_speed # already has a fallback, see above
        @eval @inline function $fun(fluid::AbstractFluid, state::NamedTuple)
            return fallback($fun, fluid, state)
        end
    end
end

# fallbacks for the state functions specific to binary fluids
for fun in setdiff(propertynames(binary_state_functions()), propertynames(all_state_functions()))
    @eval @inline function $fun(fluid::BinaryFluid, state::NamedTuple)
        return fallback($fun, fluid, state)
    end
end

# fallback implementations of potential temperature, volume, enthalpy
@inline potential_volume(fluid::SimpleFluid, state::NamedTuple)         = specific_volume(fluid, (; p=fluid.p0, s=specific_entropy(fluid, state)))
@inline potential_temperature(fluid::SimpleFluid, state::NamedTuple)    = temperature(fluid, (; p=fluid.p0, s=specific_entropy(fluid, state)))
@inline potential_enthalpy(fluid::SimpleFluid, state::NamedTuple)       = specific_enthalpy(fluid, (; p=fluid.p0, s=specific_entropy(fluid, state)))

# helper types to implement the pattern
#   fluid_pT = fluid(:p, :T)
#   volume = fluid_pT.specific_volume
#   v = @. volume(p,T)

struct FluidWithVars{vars,Fluid}
    fluid::Fluid
    function FluidWithVars{vars,Fluid}(fluid) where {vars,Fluid}
        @assert vars in state_variables(fluid)
        return new{vars,Fluid}(fluid)
    end
end

struct StateFunction{vars,fun,Fluid}
    fluid::Fluid
    function StateFunction{vars,fun,Fluid}(fluid) where {vars,fun,Fluid}
        @assert vars in state_variables(fluid)
        @assert fun in propertynames(state_functions(fluid))
        return new{vars,fun,Fluid}(fluid)
    end
end

@inlineall begin

    (fluid::SimpleFluid)(v1::Symbol, v2::Symbol) =
        FluidWithVars{(v1, v2),typeof(fluid)}(fluid)

    (fluid::BinaryFluid)(v1::Symbol, v2::Symbol, v3::Symbol) =
        FluidWithVars{(v1, v2, v3),typeof(fluid)}(fluid)

    Base.propertynames(fluid::FluidWithVars{vars}) where {vars} =
        propertynames(state_functions(fluid.fluid))

    function Base.getproperty(fluid::FluidWithVars{vars}, sym::Symbol) where {vars}
        sym in (:fluid, :vars) && return getfield(fluid, sym)
        f = getfield(fluid, :fluid)
        return StateFunction{vars,sym,typeof(f)}(f)
    end

    function (fluid_fun::StateFunction{vars,fun})(v1, v2) where {vars,fun}
        args = NamedTuple{vars}((v1, v2))
        fluid = fluid_fun.fluid
        state_fun = getfield(state_functions(fluid), fun)
        return state_fun(fluid, args)
    end
    function (fluid_fun::StateFunction{vars,fun})(v1, v2, v3) where {vars,fun}
        args = NamedTuple{vars}((v1, v2, v3))
        fluid = fluid_fun.fluid
        state_fun = getfield(state_functions(fluid), fun)
        return state_fun(fluid, args)
    end
end
