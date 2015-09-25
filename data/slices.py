from colloids.povray import *
import subprocess

for t in [120, 390, 890]:
    pos = np.loadtxt('pos_%03d'%t, skiprows=1)
    sl = (pos[:,-1]>100) & (pos[:,-1]<150)
    pos = pos[sl]
    phase = np.loadtxt('studyphase_%03d'%t, dtype=int, usecols=[1])[sl]
    
    f = File("slice_t%03d.pov"%t,"colors.inc", "header.inc")
    crystals = Union(*[
        Sphere( (x,y,z), 5)
        for x,y,z in pos[phase<3]
        ]+[Texture(Pigment(color=(0.33,0,0)))])
    f.write(crystals)
    for i, c in zip([3,4,5], [(0.66,0.33,0), (1,0.66,0.1), (0.75,0.75,1)]):
        f.write(Union(*[
            Sphere( (x,y,z), 5)
            for x,y,z in pos[phase==i]
            ]+[Texture(Pigment(color=c))
            ]))
    f.file.flush()
    subprocess.call(['povray', '-Oslice_t%03d.png'%t, '+W512', '+H512' , 'slice_t%03d.pov'%t])
