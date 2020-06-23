
struct Gene{T}
    param::Symbol
    val::T
    pow::Int
end

struct Genotype
    a::Vector{Gene}
    b::Vector{Gene}
end
