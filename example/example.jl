using Plots, BioParticleSimulation

f = BioParticleSimulation.FieldSet(100, 100)
BioParticleSimulation.fillfield(f, 0.1, 0.3, 0.00008, 0.0007) #0.0007
BioParticleSimulation.round(f)
BioParticleSimulation.initplot(f)

for i = 1:20000
    BioParticleSimulation.addresource(f, 5, 350)
    BioParticleSimulation.round(f)
end
p = BioParticleSimulation.initplot(f)


#Graph

f = BioParticleSimulation.FieldSet(100, 100)
BioParticleSimulation.fillfield(f, 0.1, 0.3, 0.00008, 0.0007)
n = Vector{Int}(undef, 0)
e = Vector{Int}(undef, 0)
anim = @animate for i ∈ 1:600
    BioParticleSimulation.addresource(f, 5, 250)
    BioParticleSimulation.round(f)
    push!(n, length(f.list))
    push!(e, f.e)
    plot(BioParticleSimulation.initplot(f), plot(n, legend = false), plot(e, legend = false),layout = grid(3, 1, heights=[0.6, 0.2, 0.2]), size=(400, 800))
end
gif(anim, "anim_fps15.gif", fps = 15)
