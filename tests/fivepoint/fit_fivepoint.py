#!/usr/bin/env python3
"""Five-point charge scan -> capacitance and grand-potential parabola, compared with the reference values below.
   fit_fivepoint.py fp_xm10.out fp_xm05.out fp_x00.out fp_x05.out fp_x10.out"""
import sys, re
HA, PHI, AREA = 27.211386245988, 4.43, 106.4615
E_TO_UC = 1.602176634e-19/(AREA*1e-16)*1e6
REF = {  # Fermi energies [Ha] of the reference run (24 ranks)
 -1.0: -0.22556788445803, -0.5: -0.20937012079592, 0.0: -0.19375886312432, 0.5: -0.17863256587943, 1.0: -0.16284938159013}
REF_FIT = dict(C=1.1764, Upzc=0.8500, R2=0.999880, Cd=17.70)
n, ef = [], []
for f in sys.argv[1:]:
    t = open(f).read()
    nex = float(re.findall(r'EXCESS_ELECTRONS\s+(\S+)', open(f.replace('.out', '.inp')).read())[0])
    e = float(re.findall(r'Fermi energy:\s+(\S+)', t)[-1])
    n.append(nex); ef.append(e)
    print(f'{f:14s} N_exc {nex:+.1f}  E_F {e:.14f}  ref {REF[nex]:.14f}  diff {e-REF[nex]:+.2e} Ha')
U = [-e*HA - PHI for e in ef]
m = len(n); sx, sy = sum(U), sum(n); sxx = sum(u*u for u in U); sxy = sum(u*v for u, v in zip(U, n))
a = (m*sxy - sx*sy)/(m*sxx - sx*sx); b = (sy - a*sx)/m
C = -a; Upzc = b/C
res = [v - (a*u + b) for u, v in zip(U, n)]
ss_res = sum(r*r for r in res); ss_tot = sum((v - sy/m)**2 for v in n); R2 = 1 - ss_res/ss_tot
print(f'fit: C = {C:.4f} e/V = {C*E_TO_UC:.2f} uF/cm^2   U_pzc = {Upzc:.4f} V vs SHE   R^2 = {R2:.6f}   max|res| = {max(abs(r) for r in res):.5f} e')
print(f'ref: C = {REF_FIT["C"]:.4f} e/V = {REF_FIT["Cd"]:.2f} uF/cm^2   U_pzc = {REF_FIT["Upzc"]:.4f} V vs SHE   R^2 = {REF_FIT["R2"]:.6f}')
print(f'parabola: dOmega(U) = -{0.5*C:.5f} (U - {Upzc:.4f})^2 eV')
worst = max(abs(e - REF[k]) for k, e in zip(n, ef))
print('VERDICT vs reference:', 'PASS - every E_F within 1e-7 Ha of the reference' if worst < 1e-7 else f'FAIL - max |dE_F| = {worst:.2e} Ha')
