RUN_FOKKER_PLANCK=0

if RUN_FOKKER_PLANCK THEN BEGIN 

;; this will run two pitch-angle models:
;; k=0 (isotropic injection), k=1 (beam field-align with a gaussian width)
for k=0,1 do begin
ss=1e9
s=interpol([-ss,ss],101)

b=bfield(s,mr=[1.1d,4d],model=1,dbds=dbds,/prec)
;; 
 
iso=([2.1,0.3])[k]
label=(['isot','beam'])[k]

;; input parameters
p=[50,99,29,200 $ ;; grid:E,mu,s,t (ned, nud, nsd, ntd)
    ,10.0,500.0,-ss,ss,2.0,4.0,2 $ ;; , emin, emax, smin, smax, mr, nexp, bmodel,
    ,1e9,1e9,5.0,3.0 $;; nmin, nmax, tmax, delta
    ,iso,2.0,1e8,0.0,1.0,0.25,1e12 $ ;; , dp, p0, dx, x0, t0, tau, lambda, 
    ,1,1,1,1,0,1] ;; flags: run_transport, run_energyloss, run_mirroring, run_diffusion,
 ;; run_turbulence, run_timeinject

;; this will call the external FP code
fp,p,r,bfield=[[s],[b]]

save,r,p,file='fp_example_'+label+'.sav'

endfor

ENDIF 

end