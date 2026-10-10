import numpy as np, time
import exp4
from run_4g import *
pos,nb,rev,d0,LB,LB0,bip=load_tiled('raw_tiled_test'); A=adj(nb); c=CT['raw_tiled_test']
print('N',len(nb),'LB',LB)
# reference: exp4 loader for the same tiled network
p2,n2,r2,d2,L2=exp4.load_net('raw_tiled'); print('same net:',np.allclose(p2%L2,pos%LB), (n2==nb).all(), (r2==rev).all(), abs(L2-LB))
A2=exp4.adj(n2)
dn='123'; u=D2[dn]; vf=0.8; g=1/np.sqrt(1-vf**2); kv=g*eta*vf/c*u
w=omega(d0,kv,eta); dk=1e-4*u; vg=(omega(d0,kv+dk,eta)-omega(d0,kv-dk,eta))/2e-4
c0=np.array((0.25,0.25,0.75))*LB0
phiA,ei,eo,info=build_packet(pos,nb,d0,LB,A,c0,kv,-vg*u,w,40.0,None)
phiB,ei2,eo2=exp4.build_packet(p2%L2,n2,d2,L2,A2,c0,kv,-vg*u,w,eta,40.0,None)
print('packet max diff',np.abs(phiA-phiB).max(),'max |phi|',np.abs(phiB).max(),'e_in/out',ei,eo,ei2,eo2)
# tracking comparison over 100 ticks, old numpy code vs compiled moments
phi=phiA.copy(); out=np.zeros_like(phi); ce=np.cos(eta)/2; em=np.exp(-1j*eta); z3=np.zeros(3); rc2=(0.45*LB)**2
cen=c0.copy(); un=c0.copy(); cn=c0.copy(); unn=c0.copy(); md=0; mp=0
for t in range(101):
    if t%10==0:
        rho=(np.abs(phi)**2).sum(1)
        for _ in range(3):
            r=pos-cen; r-=LB*np.round(r/LB); sh=(rho[:,None]*r).sum(0)/rho.sum(); cen=(cen+sh)%LB; un=un+sh
        r=pos-cen; r-=LB*np.round(r/LB); R=r+(un-c0); S=phi.sum(1); pa=np.angle((S*np.exp(-1j*(R@kv))).sum())
        for _ in range(3):
            m=moments(phi,pos,cn,LB,rc2,kv,z3); sh=m[1:4]/m[0]; cn=(cn+sh)%LB; unn=unn+sh
        m=moments(phi,pos,cn,LB,rc2,kv,unn-c0); pb=np.angle(m[9]+1j*m[10])
        md=max(md,np.abs(un-unn).max()); mp=max(mp,abs(np.angle(np.exp(1j*(pa-pb)))))
    if t<100: step(phi,nb,rev,ce,em,out); phi,out=out,phi
print('max centre diff (A)',md,'max phase diff',mp)
t1=time.time()
for _ in range(10): step(phi,nb,rev,ce,em,out); phi,out=out,phi
print('step time 1.73M int32: %.3fs'%((time.time()-t1)/10))
t1=time.time(); m=moments(phi,pos,cn,LB,rc2,kv,z3); print('moments time %.3fs'%(time.time()-t1))
