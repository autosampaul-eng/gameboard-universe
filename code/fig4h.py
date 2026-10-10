import json, glob, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=[json.load(open(f)) for f in sorted(glob.glob('res_4h/job_*.json'))]
fig,(a1,a2)=plt.subplots(1,2,figsize=(13.5,4.9))
grp=[(lambda r:r['set']=='A' and r['dir']=='100','#A84B2F','-','glass 0.8c [100]'),(lambda r:r['set']=='A' and r['dir']=='123','#20808D','-','glass 0.8c [1,2,3]'),
     (lambda r:r['set']=='D','#944454','--','glass 0.8c [100], 326 Å box'),(lambda r:r['set']=='B','#D19900','-','glass 0.5c (both directions)'),
     (lambda r:r['set']=='C','#1B474D',':','diamond (0.8c and 0.5c)')]
for f,col,ls,lab in grp:
    g=[r for r in R if f(r)]; t=np.array(g[0]['ts']); S=np.array([r['ds_series'] for r in g])
    a1.plot(t,S.mean(0),ls,color=col,lw=2,label=f'{lab} (n={len(g)})')
a1.axvspan(250,600,color='#BCE2E7',alpha=0.35,lw=0); a1.text(260,37,'scored for the scattering rate',fontsize=9,color='#28251D')
a1.set_xlabel('tick'); a1.set_ylabel('scattered share of the wave, Δs (percentage points)'); a1.set_ylim(-1,40)
a1.set_title('Glass scatters fast particles steadily from the first tick; diamond none',fontsize=10.5)
a1.legend(frameon=False,fontsize=8.5,loc='center left',bbox_to_anchor=(0.0,0.62)); a1.spines[['top','right']].set_visible(False)
for dn,col in [('100','#A84B2F'),('123','#20808D')]:
    g=[r for r in R if r['set']=='A' and r['dir']==dn]; t=np.array(g[0]['coh']['win_t'])
    a2.plot(t,np.array([r['coh']['win_v'] for r in g]).mean(0),'-',color=col,lw=2.2,label=f'{"[100]" if dn=="100" else "[1,2,3]"}: particle only')
    a2.plot(t,np.array([r['whole']['win_v'] for r in g]).mean(0),'--',color=col,lw=1.3,alpha=0.8,label=f'{"[100]" if dn=="100" else "[1,2,3]"}: whole-box centre (as in 4d–4g)')
a2.set_xlabel('tick (centre of 40-tick window)'); a2.set_ylabel('speed v/c (mean of 6 runs)'); a2.set_ylim(0.65,0.82)
a2.set_title('With the scattered waves removed, the particle moves steadily',fontsize=10.5)
a2.legend(frameon=False,fontsize=8.5,loc='upper left'); a2.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('exp4h_results.png',dpi=150)
