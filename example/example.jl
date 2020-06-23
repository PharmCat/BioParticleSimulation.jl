using Plots, BenchmarkTools, BioParticleSimulation

f = BioParticleSimulation.FieldSet(100, 100)
t = BioParticleSimulation.CellType(0, 0.3, 0.000001, 0.0007, 100, 2000, 4000, 100, 100, 10, 1)
BioParticleSimulation.addresource(f, 500, 75)
BioParticleSimulation.fillfield(f, 0.005, t) #0.0007
BioParticleSimulation.round(f)
BioParticleSimulation.plotfield(f)

for i = 1:10000

    BioParticleSimulation.round(f)
    BioParticleSimulation.growresource(f, 0.005)
    BioParticleSimulation.round(f)
end

p = BioParticleSimulation.plotfield(f)
BioParticleSimulation.round(f)
BioParticleSimulation.printall(f)


#Graph

f = BioParticleSimulation.FieldSet(100, 100)
BioParticleSimulation.fillfield(f, 0.005)
BioParticleSimulation.addresource(f, 500, 75)
n = Vector{Int}(undef, 0)
e = Vector{Float64}(undef, 0)
anim = @animate for i ∈ 1:600
    BioParticleSimulation.round(f)
    BioParticleSimulation.growresource(f, 0.005)
    BioParticleSimulation.round(f)
    push!(n, length(f.list))
    push!(e, sum(f.rfield))
    plot(BioParticleSimulation.plotfield(f), plot(n, legend = false), plot(e, legend = false),layout = grid(3, 1, heights=[0.6, 0.2, 0.2]), size=(400, 800))
end
gif(anim, "anim_fps15.gif", fps = 15)


n = Vector{Int}(undef, 0)
e = Vector{Float64}(undef, 0)
for i = 1:40000
    BioParticleSimulation.round(f)
    BioParticleSimulation.growresource(f, 0.005)
    BioParticleSimulation.round(f)
    push!(n, length(f.list))
    push!(e, sum(f.rfield))
end
plot(plot(n, legend = false), plot(e, legend = false), layout = grid(2, 1, heights=[0.5, 0.5]), size=(400, 800))
#benchmark



bench = @benchmark for i = 1:2000
    BioParticleSimulation.round(f)
    BioParticleSimulation.growresource(f, 0.005)
end
