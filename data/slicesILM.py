from colloids.povray import *
import subprocess

#4 -> LIQUID
#5 -> GAS
#3 -> SURFACE
#0 -> BCC
#1 -> HCP
#2 -> FCC

for t in [120, 390, 890]:
    pos = np.loadtxt('pos_%03d'%t, skiprows=1)
    sl = (pos[:,-1]>100) & (pos[:,-1]<150)
    pos = pos[sl]
    phase = np.loadtxt('studyphase_%03d'%t, dtype=int, usecols=[1])[sl]
    
    f = File("sliceILM_t%03d.pov"%t,"colors.inc", "header.inc")
    crystals = Union(*[
        Sphere( (x,y,z), 5)
        for x,y,z in pos[phase<3]
        ]+[Texture(Pigment(color=(0.694,  0.145,  0.623)))])
    f.write(crystals)
    for i, c in zip([3,4,5], [(0.847,0.57,0.81), (0.847,0.57,0.81), (1.,  0.455,  0.156)]):
        f.write(Union(*[
            Sphere( (x,y,z), 5)
            for x,y,z in pos[phase==i]
            ]+[Texture(Pigment(color=c))
            ]))
    f.file.flush()
    subprocess.call(['povray', '-O../presentation/sliceILM_t%03d.png'%t, '+W256', '+H256' , 'sliceILM_t%03d.pov'%t])
