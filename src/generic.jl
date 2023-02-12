for func in (:(==), :isequal, :isapprox),
    type in (:TimeEvolution, :StateSpaceRealizable, :RealizationMethod)

    @eval function Base.$func(P1::U, P2::V; kwargs...) where {U<:$type,V<:$type}
        nameof(U) === nameof(V) || return false
        fields = fieldnames(U)
        fields === fieldnames(V) || return false

        for f in fields
            isdefined(P1, f) && isdefined(P2, f) || return false
            getfield(P1, f) === getfield(P2, f) ||
                $func(getfield(P1, f), getfield(P2, f); kwargs...) ||
                return false
        end

        return true
    end
end



# think I need a IsProper / IsNotProper trait here
ssrealize(h::StateSpaceRealizable{E}) where {E} = Base.splat(ProperStateSpace{E})(ssparams(h)) 


function sqr_magnitude_response(h::StateSpaceRealizable, zs)
    n = length(zs)
    m = ninputs(h)
    mag2 = [zeros(eltype(zs), m, m) for i = 1:n]
    sqr_magnitude_response!(mag2, h, zs)
    return mag2
end

function impulse_response(h::StateSpaceRealizable, ts)
    n = length(ts)
    hs = [zeros(eltype(ts), noutputs(h), ninputs(h)) for i = 1:n]
    sqr_magnitude_response!(hs, h, ts)
    return hs
end


