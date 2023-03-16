function +(
    s1::MatrixStateSpace{E,T},
    s2::MatrixStateSpace{E,T},
    ) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    D = s2.D + s1.D
    return StateSpace{E}(A, B, C, D)
end

function +(
    s1::MatrixStateSpace{E,T},
    s2::MatrixProperStateSpace{E,T},
    ) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    D = s1.D
    return StateSpace{E}(A, B, C, D)
end

function +(
    s1::MatrixProperStateSpace{E,T},
    s2::MatrixStateSpace{E,T},
    ) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    D = s2.D
    return StateSpace{E}(A, B, C, D)
end

function +(
    s1::MatrixProperStateSpace{E,T},
    s2::MatrixProperStateSpace{E,T},
    ) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    return ProperStateSpace{E}(A, B, C)
end

-(s1::AbstractStateSpace{E,T}, s2::AbstractStateSpace{E,T}) where {E,T} = +(s1, -s2)


function ∘(
    s2::MatrixStateSpace{E,T},
    s1::MatrixStateSpace{E,T},
    ) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, s2.B * s1.D)
    C = hcat(s2.D * s1.C, s2.C)
    D = s2.D * s1.D
    return StateSpace{E}(A, B, C, D)
end


function ∘(
    s2::MatrixStateSpace{E,T},
    s1::MatrixProperStateSpace{E,T},
    ) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, zeros(T, nstates(s2), ninputs(s1)))
    C = hcat(s2.D * s1.C, s2.C)
    return ProperStateSpace{E}(A, B, C)
end


function ∘(
    s2::MatrixProperStateSpace{E,T},
    s1::MatrixStateSpace{E,T},
    ) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, s2.B * s1.D)
    C = hcat(zero(s1.C), s2.C)
    return ProperStateSpace{E}(A, B, C)
end


function ∘(
    s2::MatrixProperStateSpace{E,T},
    s1::MatrixProperStateSpace{E,T},
    ) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, zeros(T, nstates(s2), ninputs(s2)))
    C = hcat(zeros(T, noutputs(s1), nstates(s1)), s2.C)
    return ProperStateSpace{E}(A, B, C)
end