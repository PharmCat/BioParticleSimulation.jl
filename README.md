# BioParticleSimulation.jl

Simple biological particle simulation system.

Alpha version.



```
using  BioParticleSimulation

f = BioParticleSimulation.FieldSet(100, 100)
t = BioParticleSimulation.CellType(0, 0.3, 0.000001, 0.0007, 100, 2000, 4000, 100, 100, 10, 1)
BioParticleSimulation.fillfield(f, 0.1, t) #0.0007
BioParticleSimulation.round(f)
BioParticleSimulation.plotfield(f)

for i = 1:20000
    BioParticleSimulation.addresource(f, 5, 350)
    BioParticleSimulation.round(f)
    BioParticleSimulation.growresource(f, 1.005)
    BioParticleSimulation.round(f)
end

p = BioParticleSimulation.plotfield(f)
```


```

using Plots, BenchmarkTools, BioParticleSimulation

#Graph

f = BioParticleSimulation.FieldSet(100, 100)
BioParticleSimulation.fillfield(f, 0.1)
n = Vector{Int}(undef, 0)
e = Vector{Int}(undef, 0)
anim = @animate for i ∈ 1:1500
    BioParticleSimulation.addresource(f, 5, 250)
    BioParticleSimulation.round(f)
    BioParticleSimulation.growresource(f, 1.005)
    BioParticleSimulation.round(f)
    push!(n, length(f.list))
    push!(e, sum(f.rfield))
    plot(BioParticleSimulation.plotfield(f), plot(n, legend = false), plot(e, legend = false),layout = grid(3, 1, heights=[0.6, 0.2, 0.2]), size=(400, 800))
end
gif(anim, "anim_fps15.gif", fps = 15)


#benchmark

bench = @benchmark for i = 1:2000
    BioParticleSimulation.addresource(f, 5, 350)
    BioParticleSimulation.round(f)
    BioParticleSimulation.growresource(f, 1.005)
    BioParticleSimulation.round(f)
end


```
