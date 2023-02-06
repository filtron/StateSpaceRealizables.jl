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