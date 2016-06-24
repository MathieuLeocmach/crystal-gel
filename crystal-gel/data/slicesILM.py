from colloids.povray import *
import subprocess

#4 -> LIQUID
#5 -> GAS
#3 -> SURFACE
#0 -> BCC
#1 -> HCP
#2 -> FCC

for sample, times in zip(['', '_highcp'], [[120, 390, 890], [800]]):
    for t in times:
        pos = np.loadtxt('pos%s_%03d'%(sample, t), skiprows=1)
        sl = (pos[:,-1]>100) & (pos[:,-1]<125)
        pos = pos[sl]
        phase = np.loadtxt('studyphase%s_%03d'%(sample, t), dtype=int, usecols=[1])[sl]
        
        f = File("sliceILM%s_t%03d.pov"%(sample, t),"colors.inc", "header.inc")
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
        subprocess.call(['povray', '-O../presentation/sliceILM%s_t%03d.png'%(sample, t), '+W256', '+H256' , 'sliceILM%s_t%03d.pov'%(sample, t)])
