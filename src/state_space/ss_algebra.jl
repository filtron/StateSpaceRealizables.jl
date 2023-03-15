


+(s1::AbstractStateSpace{E,T}, s2::AbstractStateSpace{E,T}) where {E,T} =
ssadd(s1, ssparameter_type(s1), s2, ssparameter_type(s2))

function ssadd(
s1::StateSpace{E,T},
::MatrixParameter,
s2::StateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
B = vcat(s2.B, s1.B)
C = hcat(s2.C, s1.C)
D = s2.D + s1.D
return StateSpace{E}(A, B, C, D)
end

function ssadd(
s1::StateSpace{E,T},
::MatrixParameter,
s2::ProperStateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
B = vcat(s2.B, s1.B)
C = hcat(s2.C, s1.C)
D = s1.D
return StateSpace{E}(A, B, C, D)
end

function ssadd(
s1::ProperStateSpace{E,T},
::MatrixParameter,
s2::StateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
B = vcat(s2.B, s1.B)
C = hcat(s2.C, s1.C)
D = s2.D
return StateSpace{E}(A, B, C, D)
end

function ssadd(
s1::ProperStateSpace{E,T},
::MatrixParameter,
s2::ProperStateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
B = vcat(s2.B, s1.B)
C = hcat(s2.C, s1.C)
return ProperStateSpace{E}(A, B, C, D)
end


-(s::AbstractStateSpace) = ssnegate(s, ssparameter_type(s))
-(s1::AbstractStateSpace{E,T}, s2::AbstractStateSpace{E,T}) where {E,T} = +(s1, -s2)


ssnegate(s::ProperStateSpace{E}, ::MatrixParameter) where {E} =
ProperStateSpace{E}(s.A, s.B, -s.C)


∘(s2::AbstractStateSpace{E,T}, s1::AbstractStateSpace{E,T}) where {E,T} =
sscompose(s2, ssparameter_type(s2), s1, ssparameter_type(s1)) # prefer this one
*(s2::AbstractStateSpace{E,T}, s1::AbstractStateSpace{E,T}) where {E,T} = ∘(s2, s1)  # keep this one to not break things

function sscompose(
s2::StateSpace{E,T},
::MatrixParameter,
s1::StateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
B = vcat(s1.B, s2.B * s1.D)
C = hcat(s2.D * s1.C, s2.C)
D = s2.D * s1.D
return StateSpace{E}(A, B, C, D)
end

function sscompose(
s2::StateSpace{E,T},
::MatrixParameter,
s1::ProperStateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
B = vcat(s1.B, zeros(T, nstates(s2), ninputs(s1)))
C = hcat(s2.D * s1.C, s2.C)
return ProperStateSpace{E}(A, B, C)
end

function sscompose(
s2::ProperStateSpace{E,T},
::MatrixParameter,
s1::StateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
B = vcat(s1.B, s2.B * s1.D)
C = hcat(zero(s1.C), s2.C)
return ProperStateSpace{E}(A, B, C)
end

function sscompose(
s2::ProperStateSpace{E,T},
::MatrixParameter,
s1::ProperStateSpace{E,T},
::MatrixParameter,
) where {E,T}
A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
B = vcat(s1.B, zeros(T, nstates(s2), ninputs(s2)))
C = hcat(zeros(T, noutputs(s1), nstates(s1)), s2.C)
return ProperStateSpace{E}(A, B, C)
end