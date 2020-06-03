module BioParticleSimulation

using Random, StatsBase, Plots

export FieldSet, fillfield, initplot, addresource, round
    abstract type AbstractField end

    const UP    = ( 1,  0)
    const DOWN  = (-1,  0)
    const LEFT  = ( 0, -1)
    const RIGHT = ( 0,  1)

    mutable struct FieldSet
        field::Matrix{Union{AbstractField, Bool}}
        list::Vector{AbstractField}
        e::Int
        function FieldSet(x::Int, y::Int)::FieldSet
            new(Matrix{Union{AbstractField, Bool}}(undef, y, x), Vector{AbstractField}(undef, 0), 0)::FieldSet
        end
    end

    function Base.length(f::FieldSet)
        length(f.list)
    end

    mutable struct Cell <: AbstractField
        enrg::Int
        x::Int
        y::Int
        movep::Float32
        killp::Float32
        replp::Float32
        field::FieldSet
        maxe::Int
        maxl::Int
        function Cell(enrg, x, y, movep, killp, replp, field, maxe)::Cell
            new(enrg, x, y, movep, killp, replp, field, maxe, 4000)::Cell
        end
        function Cell(enrg, x, y, movep, killp, replp, field, maxe, maxl)::Cell
            new(enrg, x, y, movep, killp, replp, field, maxe, maxl)::Cell
        end
    end


    struct Resource <: AbstractField
        val
    end


    function fillfield(f, p, m, k, rep)
        x    = size(f.field, 2)
        y    = size(f.field, 1)
        mask = sample([true, false], Weights([p, 1-p]), (x,y))
        for r = 1:y
            for c = 1:x
                if mask[r, c]
                    newcell = Cell(1000, c, r, m, k, rep, f, 2000)
                    push!(f.list, newcell)
                    f.field[r, c] = newcell
                else
                    f.field[r, c] = false
                end
            end
        end
    end

    function addresource(f, n, val)
        x     = size(f.field, 2)
        y     = size(f.field, 1)

        for i = 1:n
            rx = rand(1:x)
            ry = rand(1:y)
            if emptycellonfield(f, rx, ry)
                f.field[ry, rx] = Resource(val)
                f.e += 1
            end
        end
    end

    function cellonfield(f, x, y)
        isa(f.field[y, x], Cell)
    end
    function resonfield(f, x, y)
        isa(f.field[y, x], Resource)
    end
    function emptycellonfield(f, x, y)
        if isa(f.field[y, x], Bool)
            return true
        elseif isa(f.field[y, x], Cell)
            if f.field[y, x].enrg <= 0 || f.field[y, x].maxl <= 0 return true end
        end
        false
    end

    function movecell(c)
        vec = sample([UP, DOWN, LEFT, RIGHT])
        newx = c.x + vec[1]
        newy = c.y + vec[2]
        if c.enrg > 10 && c.maxl > 0 && newx > 0 && newy > 0 && newx <= size(c.field.field, 2) && newy <= size(c.field.field, 1)
            if cellonfield(c.field, newx, newy)
                killcell(c, newx, newy)
                return
            end
            c.enrg -= 10
            if resonfield(c.field, newx, newy)
                if c.enrg + getfield(c, newx, newy).val >= c.maxe
                    c.enrg = c.maxe
                else
                    c.enrg += getfield(c, newx, newy).val
                end
                c.field.e -= 1
            end
            setcellfield(c, false)
            c.x = newx
            c.y = newy
            setcellfield(c, c)
        else
            idlecell(c)
        end
    end
    function replcell(c)
        vec = sample([UP, DOWN, LEFT, RIGHT])
        newx = c.x + vec[1]
        newy = c.y + vec[2]
        if c.enrg > 100 && c.maxl > 0 && newx > 0 && newy > 0 && newx <= size(c.field.field, 2) && newy <= size(c.field.field, 1) && emptycellonfield(c.field, newx, newy)
            c.enrg -= 100
            newcell = Cell(100, newx, newy, c.movep, c.killp, c.replp, c.field, 2000)
            push!(c.field.list, newcell)
            c.field.field[newy, newx] = newcell
        else
            idlecell(c)
        end
    end

    function killcell(c, x, y)
        if sample([true, false], Weights([c.killp, 1 - c.killp]))
            setcellfield(c, false)
            if c.enrg + 100 >= c.maxe
                c.enrg = c.maxe
            else
                c.enrg += 100
            end
            getfield(c, x, y).enrg = 0
            getfield(c, x, y).maxl = 0
            c.x = x
            c.y = y
            setcellfield(c, c)
        else
            idlecell(c)
        end
    end
    function idlecell(c)
        c.enrg -= 1
        if c.enrg <= 0
            setcellfield(c, false)
        end
    end

    function setcellfield(c, val)
        c.field.field[c.y, c.x] = val
    end
    function setfield(c, x, y, val)
        c.field.field[y, x] = val
    end
    function getfield(c, x, y)
        c.field.field[y, x]
    end

    function cellaction(c::Cell)
        if c.enrg > 0 && c.maxl > 0
            action =  sample([movecell, replcell, idlecell], Weights([c.movep, c.replp, 1 - c.movep - c.replp]))
            action(c)
            c.maxl -= 1
            if c.maxl <= 0
                setcellfield(c, false)
            end
        end
    end

    function cellaction(r::Resource)
    end
    function cellaction(b::Bool)
    end


    function round(f::FieldSet)
        #=
        x    = size(f.field, 2)
        y    = size(f.field, 1)
        for r = 1:y
            for c = 1:x
                cellaction(f.field[r,c])
            end
        end
        =#
        for c in f.list
            cellaction(c)
        end
        deleteat!(f.list, findall(x -> x.enrg <=0 || x.maxl <= 0, f.list))
        nothing
    end

    function round(f::FieldSet, n::Int)
        for i = 1:n
            round(f)
        end
    end

    function initplot(f)
        cnt  = 0
        cntr = 0
        #=
        for r = 1:size(f.field, 1)
            for c = 1:size(f.field, 2)
                if isa(f.field[r, c], Cell) cnt += 1 end
                if isa(f.field[r, c], Resource) cntr += 1 end
            end
        end
        =#
        mx   = Matrix{Int}(undef, length(f.list), 2)
        mxr  = Matrix{Int}(undef, f.e, 2)
        ccnt  = 1
        ccntr = 1
        for r = 1:size(f.field, 1)
            for c = 1:size(f.field, 2)
                if isa(f.field[r, c], Cell)
                    mx[ccnt, 1] = r
                    mx[ccnt, 2] = c
                    ccnt += 1
                end
                if isa(f.field[r, c], Resource)
                    mxr[ccntr, 1] = r
                    mxr[ccntr, 2] = c
                    ccntr += 1
                end
            end
        end
        p = scatter(mx[:,1], mx[:,2], title = "Cells", color = :red, legend = false)
        p = scatter!(p, mxr[:,1], mxr[:,2], title = "Cells", color = :green, legend = false)
        p
    end

    function Base.show(io::IO, f::FieldSet)
        x    = size(f.field, 2)
        y    = size(f.field, 1)
        for r = 1:y
            for c = 1:x
                if f.field[r, c] == false
                    print(io, " * ")
                elseif isa(f.field[r, c], Cell)
                    print(io, " C ")
                elseif isa(f.field[r, c], Resource)
                    print(io, " R ")
                else
                    print(io, "   ")
                end
            end
            println(io, "")
        end
    end

    function printall(f)
        x    = size(f.field, 2)
        y    = size(f.field, 1)
        cnt  = 0
        for r = 1:y
            for c = 1:x
                if isa(f.field[r, c], Cell)
                    println("X: $c ($(f.field[r, c].x)) Y: $r ($(f.field[r, c].y)) Energy: $(f.field[r, c].enrg) Live: $(f.field[r, c].maxl)")
                    cnt += 1
                end
            end
        end
        println("Total: $cnt")
    end

end # module
