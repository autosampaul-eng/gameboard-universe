import json, glob, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S={'s3 (Exp 4i)':(sorted(glob.glob('res_4i/job_*.json')),json.load(open('fit_4i.json')),'#20808D','o'),
   's4 (replication)':(sorted(glob.glob('res_4i_rep/job_*.json')),json.load(open('fit_4i_rep.json')),'#A84B2F','s')}
fig,(a1,a2)=plt.subplots(1,2,figsize=(13.5,5.0))
for lab,(files,F,col,m) in S.items():
    R=[json.load(open(f)) for f in files]; k=np.array(F['k0']); y=np.array(F['inv_ell']); e=np.array(F['se']); lam=2*np.pi/k/2.735
    off=0.985 if m=='o' else 1.015
    a1.errorbar(lam*off,1/y/2.735,yerr=(e/y)/y/2.735,fmt=m,color=col,ms=7,capsize=3,label=f"{lab}: exponent {F['p']:.2f} ± {F['p_se']:.2f}")
    kk=np.linspace(k.min()*0.95,k.max()*1.05,50); a1.plot(2*np.pi/kk/2.735,1/np.exp(F['lnA']+F['p']*np.log(kk))/2.735,'-',color=col,lw=1.2,alpha=0.8)
    rat=[];er=[]
    for p in range(1,6):
        g=[r for r in R if r['set']=='G' and r['point']==p]
        A=np.array([r['inv_ell'] for r in g if r['dir']=='100']); B=np.array([r['inv_ell'] for r in g if r['dir']=='123'])
        q=A.mean()/B.mean(); rat.append(100*(q-1)); er.append(100*q*np.hypot(A.std(ddof=1)/np.sqrt(3)/A.mean(),B.std(ddof=1)/np.sqrt(3)/B.mean()))
    a2.errorbar(lam*off,rat,yerr=er,fmt=m+'-',color=col,ms=7,capsize=3,lw=1.4,label=lab)
a1.set_xscale('log'); a1.set_yscale('log'); a1.set_xticks([18,21,25,30,36]); a1.set_xticklabels(['18','21','25','30','36']); a1.minorticks_off()
a1.set_yticks([500,1000,2000,4000,8000]); a1.set_yticklabels(['500','1,000','2,000','4,000','8,000'])
a1.set_xlabel('wavelength (bonds)'); a1.set_ylabel('mean free path (bonds)'); a1.legend(frameon=False,fontsize=9,loc='upper left')
a1.set_title('Two independent glasses give the same scattering law',fontsize=10.5); a1.spines[['top','right']].set_visible(False)
a2.axhline(0,color='#7A7974',lw=0.8); a2.set_xscale('log'); a2.set_xticks([18,21,25,30,36]); a2.set_xticklabels(['18','21','25','30','36']); a2.minorticks_off()
a2.set_xlabel('wavelength (bonds)'); a2.set_ylabel('scattering along [100] vs [1,2,3] (% difference)')
a2.set_title('The direction difference changes sign between samples',fontsize=10.5); a2.legend(frameon=False,fontsize=9,loc='upper left'); a2.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('exp4i_rep_results.png',dpi=150)
