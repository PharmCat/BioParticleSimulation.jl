# BioParticleSimulation.jl

Simple biological particle simulation system.

Alpha version.

```
using Plots, BioParticleSimulation

f = BioParticleSimulation.FieldSet(100, 100)
BioParticleSimulation.fillfield(f, 0.1, 0.3, 0.00008, 0.0007)
BioParticleSimulation.round(f)
BioParticleSimulation.initplot(f)

for i = 1:20000
    BioParticleSimulation.addresource(f, 5, 350)
    BioParticleSimulation.round(f)
end
p = BioParticleSimulation.initplot(f)
```
