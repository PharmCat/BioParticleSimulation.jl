module BioParticleSimulation

using Random, StatsBase, Plots

export FieldSet, fillfield, plotfield, addresource, round
    abstract type AbstractField end

    const UP    = ( 1,  0)
    const DOWN  = (-1,  0)
    const LEFT  = ( 0, -1)
    const RIGHT = ( 0,  1)

    mutable struct FieldSet
        field::Matrix{Union{AbstractField, Bool}}
        list::Vector{AbstractField}
        rfield::Matrix{Int}
        makecell::Function
        function FieldSet(x::Int, y::Int)::FieldSet
            new(Matrix{Union{AbstractField, Bool}}(undef, y, x), Vector{AbstractField}(undef, 0), zeros(Int, y, x), makecell_simple)::FieldSet
        end
    end

    function Base.length(f::FieldSet)
        length(f.list)
    end

struct CellType
    id::Int
    movep::Float32
    killp::Float32
    replp::Float32
    kille::Int
    maxe::Int
    maxl::Int
    replcost::Int
    replinite::Int
    movecost::Int
    idlecost::Int
    function CellType(
        id,
        movep,
        killp,
        replp,
        kille,
        maxe,
        maxl,
        replcost,
        replinite,
        movecost,
        idlecost,
    )::CellType
        new(
            id,
            movep,
            killp,
            replp,
            kille,
            maxe,
            maxl,
            replcost,
            replinite,
            movecost,
            idlecost,
        )::CellType
    end
    function CellType(;
        id = 0,
        movep = 0.3,
        killp = 0.000001,
        replp = 0.0007,
        kille = 100,
        maxe = 2000,
        maxl = 4000,
        replcost = 100,
        replinite = 100,
        movecost = 10,
        idlecost = 1,
    )::CellType
        new(
            id,
            movep,
            killp,
            replp,
            kille,
            maxe,
            maxl,
            replcost,
            replinite,
            movecost,
            idlecost,
        )::CellType
    end
end

    const CELL_A = CellType(0, 0.3, 0.000001, 0.0007, 100, 2000, 4000, 100, 100, 10, 1)

    mutable struct Cell <: AbstractField
        enrg::Int
        x::Int
        y::Int
        livetime::Int
        field::FieldSet
        type::CellType
        function Cell(enrg, x, y, l, field)::Cell
            new(enrg, x, y, l, field, CELL_A)::Cell
        end
        function Cell(enrg, x, y, l, field, type)::Cell
            new(enrg, x, y, l, field, type)::Cell
        end
    end

    function fillfield(f::FieldSet, p, type)
        x    = size(f.field, 2)
        y    = size(f.field, 1)
        mask = sample([true, false], Weights([p, 1-p]), (x,y))
        for r = 1:y
            for c = 1:x
                if mask[r, c]
                    newcell = Cell(1000, c, r, type.maxl, f, type)
                    push!(f.list, newcell)
                    f.field[r, c] = newcell
                else
                    f.field[r, c] = false
                end
            end
        end
    end
    function fillfield(f::FieldSet, p)
        fillfield(f::FieldSet, p, CELL_A)
    end
"""
    fillfield(f::FieldSet)

p - cell contamination probability;
type -
"""
    function fillfield(f::FieldSet; p, type)
        fillfield(f::FieldSet, p, type)
    end

    function addresource(f::FieldSet, n, val)
        x     = size(f.rfield, 2)
        y     = size(f.rfield, 1)
        for i = 1:n
            rx = rand(1:x)
            ry = rand(1:y)
            if f.rfield[ry, rx] < 2000 f.rfield[ry, rx] += val end
            #=
            if emptycellonfield(f, rx, ry)
                f.field[ry, rx] = Resource(val)
                f.e += 1
            end
            =#
        end
    end
    function growresource(f::FieldSet, m)
        for r = 1:size(f.rfield, 1)
            for c = 1:size(f.rfield, 2)
                f.rfield[r, c] = Int(ceil(f.rfield[r, c] * m))
                if f.rfield[r, c] > 2000 f.rfield[r, c] = 2000 end
            end
        end
    end

    function cellonfield(f, x, y)
        if isa(f.field[y, x], Cell)
            if f.field[y, x].enrg > 0 && f.field[y, x].livetime > 0 return true end
        end
        false
    end
    #=
    function resonfield(f, x, y)
        isa(f.field[y, x], Resource)
    end
    =#
    function emptycellonfield(f, x, y)
        if isa(f.field[y, x], Bool)
            return true
        elseif isa(f.field[y, x], Cell)
            if f.field[y, x].enrg <= 0 || f.field[y, x].livetime <= 0 return true end
        end
        false
    end

    function movecell(c)
        vec = sample([UP, DOWN, LEFT, RIGHT])
        newx = c.x + vec[1]
        newy = c.y + vec[2]
        if c.enrg > c.type.movecost && c.livetime > 0 && newx > 0 && newy > 0 && newx <= size(c.field.field, 2) && newy <= size(c.field.field, 1)
            if cellonfield(c.field, newx, newy)
                killcell(c, newx, newy)
                return
            end
            c.enrg -= c.type.movecost

            if c.enrg + getrfield(c, newx, newy) >= c.type.maxe
                c.enrg = c.type.maxe
            else
                c.enrg += getrfield(c, newx, newy)
            end

            setrfield(c, newx, newy, 0)

            setcellfield(c, false)
            c.x = newx
            c.y = newy
            setcellfield(c, c)
        else
            idlecell(c)
        end
    end
    function replcell(c)
        if c.enrg > c.type.replcost && c.livetime > 0
            c.field.makecell(c)
        else
            idlecell(c)
        end
    end
    function makecell_simple(cp::Cell)
        vec = sample([UP, DOWN, LEFT, RIGHT])
        newx = cp.x + vec[1]
        newy = cp.y + vec[2]
        if newx > 0 && newy > 0 && newx <= size(cp.field.field, 2) && newy <= size(cp.field.field, 1) && emptycellonfield(cp.field, newx, newy)
            cp.enrg -= cp.type.replcost
            newcell = Cell(cp.type.replinite, newx, newy, cp.type.maxl, cp.field, cp.type)
            push!(cp.field.list, newcell)
            cp.field.field[newy, newx] = newcell
        else
            idlecell(cp)
        end
    end


    function killcell(c, x, y)
        if sample([true, false], Weights([c.type.killp, 1 - c.type.killp]))
            setcellfield(c, false)
            if c.enrg + getfield(c, x, y).type.kille >= c.type.maxe
                c.enrg = c.type.maxe
            else
                c.enrg += getfield(c, x, y).type.kille
            end
            getfield(c, x, y).enrg = 0
            getfield(c, x, y).livetime = 0
            c.x = x
            c.y = y
            setcellfield(c, c)
        else
            idlecell(c)
        end
    end
    function idlecell(c)
        c.enrg -= c.type.idlecost
        if c.enrg <= 0
            setcellfield(c, false)
        end
    end

    function setcellfield(c::Cell, val)
        c.field.field[c.y, c.x] = val
    end
    function setfield(c::Cell, x, y, val)
        c.field.field[y, x] = val
    end
    function getfield(c::Cell, x, y)
        c.field.field[y, x]
    end
    function getrfield(c::Cell, x, y)
        c.field.rfield[y, x]
    end
    function setrfield(c::Cell, x, y, val)
        c.field.rfield[y, x] = val
    end


    function cellaction(c::Cell)
        if c.enrg > 0 && c.livetime > 0
            action =  sample([movecell, replcell, idlecell], Weights([c.type.movep, c.type.replp, 1 - c.type.movep - c.type.replp]))
            action(c)
            c.livetime -= 1
            if c.livetime <= 0
                setcellfield(c, false)
            end
        end
    end

    #function cellaction(b::Bool)
    #end


    function round(f::FieldSet)
        for c in f.list
            cellaction(c)
        end
        deleteat!(f.list, findall(x -> x.enrg <=0 || x.livetime <= 0, f.list))
        nothing
    end

    function round(f::FieldSet, n::Int)
        for i = 1:n
            round(f)
        end
    end

    function plotfield(f)
        cnt  = 0
        mx   = Matrix{Int}(undef, length(f.list), 2)
        ccnt  = 1

        for r = 1:size(f.field, 1)
            for c = 1:size(f.field, 2)
                if isa(f.field[r, c], Cell)
                    mx[ccnt, 1] = r
                    mx[ccnt, 2] = c
                    ccnt += 1
                end
            end
        end
        p = heatmap(f.rfield, c = :greens, legend = false)
        p = scatter!(mx[:,2], mx[:,1], title = "Cells", color = :red, legend = false)
        #p = scatter!(p, mxr[:,1], mxr[:,2], title = "Cells", color = :green, legend = false)
        p
    end
    #=
    function Base.show(io::IO, f::FieldSet)
        x    = size(f.field, 2)
        y    = size(f.field, 1)
        for r = 1:y
            for c = 1:x
                if f.field[r, c] == false
                    print(io, " * ")
                elseif isa(f.field[r, c], Cell)
                    print(io, " C ")
                else
                    print(io, "   ")
                end
            end
            println(io, "")
        end
    end
    =#
    function printall(f)
        x    = size(f.field, 2)
        y    = size(f.field, 1)
        cnt  = 0
        for r = 1:y
            for c = 1:x
                if isa(f.field[r, c], Cell)
                    println("X: $c ($(f.field[r, c].x)) Y: $r ($(f.field[r, c].y)) Energy: $(f.field[r, c].enrg) Live: $(f.field[r, c].livetime)")
                    cnt += 1
                end
            end
        end
        println("Total: $cnt")
    end

end # module
