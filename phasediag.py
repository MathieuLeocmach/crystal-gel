from matplotlib.pyplot import *
from colloids.phase import *
from os.path import isfile

qR = 0.071
cpov = 1
q = qR2q(qR)

#experimental points
phixp = np.array([0.25]*4+[0.1]*2)
cpxp = np.array([0.38, 0.48, 0.57, 1.07, 0.82, 1.36])

#load or compute theoretical phase diagram
fc, pivc = Liu().critical_point(q)
if isfile('phasediag.txt'):
    GL = np.loadtxt('phasediag.txt')
else:
    GL = Liu().all_GL(q, maxpiv=3500)
    np.savetxt('phasediag.txt', GL)
if isfile('phasediagFS.txt'):
    FS = np.loadtxt('phasediagFS.txt')
else:
    FS = all_coexistence(q, Liu(), Hall(), maxpiv=3500)
    np.savetxt('phasediagFS.txt', FS)


fig = figure('experimental phase diagram')
clf()

#draw experimental points
scatter(phixp, cpxp, marker='s')

#draw critical point
scatter(f2vf(fc), piv2y(pivc, q) * alpha(fc, qR) * cpov, c='r', marker='o')

#draw binodal and spinodals
for f,l in zip(GL.T[1:], ['--']*2+['-']*2):
    cp = piv2y(GL[:,0], q) * alpha(f, qR)
    plot(f2vf(f), cp * cpov, 'k'+l)
#draw tie lines
for piv, fG, fL in GL[[50,60,70,80,85,88]+range(90,98),:3]:
    cpG = piv2y(piv, q) * alpha(fG, qR) * cpov
    cpL = piv2y(piv, q) * alpha(fL, qR) * cpov
    plot(f2vf(np.array([fG,fL])), [cpG, cpL], '-', color=(0.5,0.5,0.5))

#draw equicomposition lines
for x in np.arange(1,10)/10.:
    phiG, phiL = f2vf(GL[:,1:3]).T
    phi = x*phiG+(1-x)*phiL
    cp = piv2y(GL[:,0], q) * alpha(vf2f(phi), qR)
    plot(phi, cp * cpov, '-', color=(0.5,0.5,0.5))
    
#draw fluid-crystal
for f in FS.T[1:]:
    cp = piv2y(FS[:,0], q) * alpha(f, qR)
    plot(f2vf(f), cp * cpov, 'k.')

ylabel(r'$c_p$ (g/L)')
xlabel(r'$\phi$ (%)')
ylim(0,2)
xlim(0,0.75)
draw()
savefig('phasediag.pdf')

#save binodal and spinodal
for f,l in zip(GL.T[1:], ['bg', 'bl', 'sg', 'sl']):
    cp = piv2y(GL[:,0], q) * alpha(f, qR)
    np.savetxt(
        'gasliquid_%s.phd'%l,
        np.column_stack((f2vf(f), cp * cpov)),
        header='phi\tcp',
        )
#save fluid-crystal
for f,l in zip(FS.T[1:], ['f', 'x']):
    cp = piv2y(FS[:,0], q) * alpha(f, qR)
    np.savetxt(
        'fluidcrystal_%s.phd'%l,
        np.column_stack((f2vf(f), cp * cpov)),
        header='phi\tcp',
        )

fig = figure('theoretical phase diagram')
clf()

def xp2th(cp, phi, qR=0.1, cpov=1.):
    """Converts experimental polymer concentration in osmotic pressure"""
    return y2piv(cp/cpov/alpha(vf2f(phi), qR2q(qR)), qR)

#draw experimental points
scatter(phixp, xp2th(cpxp, phixp, qR, cpov))

#draw binodal and spinodals
for f,l in zip(GL.T[1:], ['--']*2+['-']*2):
    plot(f2vf(f), GL[:,0], 'k'+l)

#draw equicomposition lines
for x in np.arange(1,10)/10.:
    phiG, phiL = f2vf(GL[:,1:3]).T
    plot(x*phiG+(1-x)*phiL, GL[:,0], '-', color=(0.5,0.5,0.5))

#for sb, l in zip('sb', ['-','--']):
 #   for ph in 'gl':
  #      phi, cp = np.loadtxt('gasliquid_%s%s.phd'%(sb,ph), unpack=True)
   #     plot(phi*100, xp2th(cp, phi, 0.1), l+'k')
ylabel(r'$\Pi v$')
xlabel(r'$\phi$ (%)')
ylim(0,3500)
xlim(0,0.74)
draw()
savefig('phasediagth.pdf')

