from colloids.povray import *
import subprocess, os
from matplotlib.pyplot import imread, imsave

#4 -> LIQUID
#5 -> GAS
#3 -> SURFACE
#0 -> BCC
#1 -> HCP
#2 -> FCC
label = '_article'

for sample, times in zip(['', '_highcp'], [[120, 390, 890], [800]]):
    for t in times:
        pos = np.loadtxt('pos%s_%03d'%(sample, t), skiprows=1)
        sl = (pos[:,-1]>100) & (pos[:,-1]<125)
        pos = pos[sl]
        phase = np.loadtxt('studyphase%s_%03d'%(sample, t), dtype=int, usecols=[1])[sl]
        
        f = File("slice%s%s_t%03d.pov"%(label, sample, t),"colors.inc", "header.inc")
        crystals = Union(*[
            Sphere( (x,y,z), 5)
            for x,y,z in pos[phase<3]
            ]+[Texture(Pigment(color=(0.694,  0.145,  0.623)))])
        f.write(crystals)
        for i, c in zip([3,4,5], [(0.347,0.07,0.311), (0.5,0.5,0.5), (1.,  0.455,  0.156)]):
            f.write(Union(*[
                Sphere( (x,y,z), 5)
                for x,y,z in pos[phase==i]
                ]+[Texture(Pigment(color=c))
                ]))
        f.file.flush()
        jpgname = '../presentation/slice%s%s_t%03d.jpg'%(label, sample, t)
        pngname = '../presentation/slice%s%s_t%03d.png'%(label, sample, t)
        subprocess.call(['povray', '-O%s'%pngname, '+W256', '+H256' , 'slice%s%s_t%03d.pov'%(label, sample, t)])
        imsave(jpgname, imread(pngname))
        os.remove(pngname)
        
        
