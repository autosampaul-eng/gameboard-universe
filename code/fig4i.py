import json, glob, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=[json.load(open(f)) for f in sorted(glob.glob('res_4i/job_*.json'))]; F=json.load(open('fit_4i.json'))
fig,(a1,a2)=plt.subplots(1,2,figsize=(13.5,5.0))
cols=['#A84B2F','#944454','#D19900','#20808D','#1B474D']
for p in range(1,6):
    g=[r for r in R if r['set']=='G' and r['point']==p]; t=np.array(g[0]['ts'])
    S=np.array([r['ds_series'] for r in g])
    lab=f"{g[0]['k0'] and 2*np.pi/g[0]['k0']/2.735:.1f} bonds (launch {g[0]['vf']:.2f}c)"
    a1.plot(t,S.mean(0),color=cols[p-1],lw=2,label=lab)
c=[r for r in R if r['set']=='C']; a1.plot(t,np.mean([r['ds_series'] for r in c],0),':',color='#28251D',lw=1.5,label='diamond, all 5 (none)')
a1.set_xlabel('tick'); a1.set_ylabel('scattered share of the wave, Δs (percentage points)'); a1.set_ylim(-1,40)
a1.set_title('Longer wavelengths scatter far less (mean of 6 glass runs each)',fontsize=10.5); a1.legend(frameon=False,fontsize=8.5,title='wavelength',title_fontsize=8.5,loc='upper left')
a1.spines[['top','right']].set_visible(False)
k=np.array(F['k0']); y=np.array(F['inv_ell']); e=np.array(F['se']); lam=2*np.pi/k/2.735
for dn,m,off in [('100','o',-0.012),('123','s',0.012)]:
    for p in range(1,6):
        v=[r['inv_ell'] for r in R if r['set']=='G' and r['point']==p and r['dir']==dn]
        a2.plot(lam[p-1]*(1+off)*np.ones(len(v)),1/np.array(v)/2.735,m,ms=4,color='#7A7974',alpha=0.6,mfc='none',label=None)
a2.errorbar(lam,1/y/2.735,yerr=(e/y)*(1/y)/2.735,fmt='o',color='#20808D',ms=7,capsize=3,label='mean of 6 runs (± SE)')
kk=np.linspace(k.min()*0.95,k.max()*1.05,50); a2.plot(2*np.pi/kk/2.735,1/np.exp(F['lnA']+F['p']*np.log(kk))/2.735,'-',color='#A84B2F',lw=1.4,label=f"fit: mean free path ∝ wavelength^{F['p']:.2f} (±{F['p_se']:.2f})")
a2.plot([],[],'o',color='#7A7974',mfc='none',ms=4,label='single runs: ○ [100], □ [1,2,3]')
a2.set_xscale('log'); a2.set_yscale('log'); a2.set_xlabel('wavelength (bonds)'); a2.set_ylabel('mean free path (bonds)')
a2.set_xticks([18,21,25,30,36]); a2.set_xticklabels(['18','21','25','30','36']); a2.minorticks_off()
a2.set_yticks([500,1000,2000,4000,8000]); a2.set_yticklabels(['500','1,000','2,000','4,000','8,000'])
a2.set_title('Doubling the wavelength makes the mean free path 13× longer',fontsize=10.5); a2.legend(frameon=False,fontsize=8.5,loc='upper left')
a2.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('exp4i_results.png',dpi=150)
print(1/y[0]/2.735, 1/y[-1]/2.735, y[0]/y[-1])
