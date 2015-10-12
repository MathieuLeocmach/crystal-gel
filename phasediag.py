from matplotlib.pyplot import *
from colloids.phase import *
from os.path import isfile
from scipy.constants import N_A

#radius of gyration in nm
R = 82.
#Molecular weight (Dalton or g/mol)
Mw = 3.6e6
#colloid diameter in nm
sigma = 2200.
#density of the solvent in g/mL
d = 1.24

qR = 2*R/sigma
cpov = 1000 * (3 * Mw)/(4 * np.pi * (R*1e-9)**3 * N_A * d*1e6)
q = qR2q(qR)
print 'qR = %0.3f, q = %0.3f, Cov = %0.3f mg/g'%(qR, q, cpov)

#experimental points
phixpG = np.array([0.069, 0.128, 0.135, 0.126, 0.333, 0.313, 0.294, 0.305])
cpxpG = np.array([0.9, 0.786, 0.82, 1.36, 1.07, 0.38, 0.57, 0.48])
phixpF = np.array([0.141, 0.141, 0.139, 0.084, 0.307, 0.309])
cpxpF = np.array([0.25,0.48,0.57, 0.76, 0.28, 0.34])

#load or compute theoretical phase diagram
fc, pivc = Liu().critical_point(q)
print 'At critical point colloid volume fraction is %0.3f and osmotic insertion work is %0.3f kT'%(f2vf(fc), pivc)
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
scatter(phixpG, cpxpG, marker='s')
scatter(phixpF, cpxpF, color='k')

#draw critical point
scatter(f2vf(fc), piv2y(pivc, qR) * alpha(fc, q) * cpov, c='r', marker='o')

#draw and save binodal and spinodals
for f,l,ph in zip(GL.T[1:], ['--']*2+['-']*2, ['bg', 'bl', 'sg', 'sl']):
    y = piv2y(GL[:,0], qR) * alpha(f, q)
    plot(f2vf(f), y * cpov, 'k'+l)
    np.savetxt(
        'gasliquid_%s.phd'%ph,
        np.column_stack((f2vf(f), y * cpov)),
        header='phi\tcp',
        )
#draw tie lines
for piv, fG, fL in GL[[50,60,70,80,85,88]+range(90,98),:3]:
    cpG = piv2y(piv, qR) * alpha(fG, q) * cpov
    cpL = piv2y(piv, qR) * alpha(fL, q) * cpov
    plot(f2vf(np.array([fG,fL])), [cpG, cpL], '-', color=(0.5,0.5,0.5))

#draw equicomposition lines
for x in np.arange(1,10)/10.:
    phiG, phiL = f2vf(GL[:,1:3]).T
    phi = x*phiG+(1-x)*phiL
    cp = piv2y(GL[:,0], qR) * alpha(vf2f(phi), q)
    plot(phi, cp * cpov, '-', color=(0.5,0.5,0.5))
    
#draw and save fluid-crystal
for f, ph in zip(FS.T[1:], 'fx'):
    y = piv2y(FS[:,0], qR) * alpha(f, q)
    plot(f2vf(f), y * cpov, 'k.')
    np.savetxt(
        'fluidcrystal_%s.phd'%ph,
        np.column_stack((f2vf(f), y * cpov)),
        header='phi\tcp',
        )

ylabel(r'$c_p$ (mg/g)')
xlabel(r'$\phi$ (%)')
ylim(0,2)
xlim(0,0.75)
draw()
savefig('phasediag.pdf')

fig = figure('theoretical phase diagram')
clf()

def xp2th(cp, phi, qR=0.1, cpov=1.):
    """Converts experimental polymer concentration in osmotic pressure"""
    return y2piv(cp/cpov/alpha(vf2f(phi), qR2q(qR)), qR)

#draw experimental points
scatter(phixpG, xp2th(cpxpG, phixpG, qR, cpov), marker='s')
scatter(phixpF, xp2th(cpxpF, phixpF, qR, cpov), c='k')

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

