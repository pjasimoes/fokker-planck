FUNCTION bfield,si,mr=mr,n=n,b0=b0,model=model,dbds=dbds,precipitation=precipitation

 if n_elements(mr) eq 0 then mr=1.0d
 if n_elements(n) eq 0 then n=2.0d
 if n_elements(b0) eq 0 then b0=1.0d
 if n_elements(model) eq 0 then model=1

if n_elements(mr) gt 1 then BEGIN
  b1=bfield(si,mr=mr[0],model=model,dbds=dbds1,precipitation=precipitation)
  b2=bfield(si,mr=mr[1],model=model,dbds=dbds2,precipitation=precipitation)
  b=[b1[0:n_elements(si)/2-1],b2[n_elements(si)/2:*]]
  dbds=[dbds1[0:n_elements(si)/2-1],dbds2[n_elements(si)/2:*]]
  return,b

endif 

 if ~keyword_set(precipitation) then BEGIN 

    l=si[n_elements(si)-1] 

    CASE (model) OF 

       1: BEGIN 
          ;; Parabolic model: 
          ;; B(s) = B0 * (1+s^n/LB^n)  
          ;; where: LB = l / sqrt(mr-1.0)
          ;; l: half loop length (cm)
          ;; B0: magnetic field at looptop B(l/2)=B0
          ;; mr: mirror ratio B(s)/B0
          ;; n: n=2 (parabolic), higher value can be used (n=4)
          print,'parabolic model'
          LB = l / sqrt(mr-1.0)                         
          b = b0 * (1.0 + si^n / LB^n)
          dbds = n * si^(n-1.0) / ((LB^n) + (si^n)) 
       END

       2: BEGIN  
          ;; Exponential model: 
          ;; B(s) = B0 * exp(log(mr)*(s/l)^n)
          ;; where:
          ;; L: half loop length (cm)
          ;; B0: magnetic field at looptop B(L/2)=B0
          ;; mr: mirror ratio B(s)/B0
          ;; n: higher values cause stronger convergence near the footpoints
          print,'exponential model'
          dbds = alog(mr) / (l^n) * n * si^(n-1.0) ; // d (ln B) / ds
          b = b0 * exp(alog(mr) * (si/l)^n)        ;           // b(s)
       END 
    ENDCASE
    
    return, b

 ENDIF ELSE BEGIN 

    l=si[n_elements(si)-2] 

    CASE (model) OF 

       1: BEGIN 
          ;; Parabolic model: 
          ;; B(s) = B0 * (1+s^n/LB^n)  
          ;; where: LB = l / sqrt(mr-1.0)
          ;; l: half loop length (cm)
          ;; B0: magnetic field at looptop B(l/2)=B0
          ;; mr: mirror ratio B(s)/B0
          ;; n: n=2 (parabolic), higher value can be used (n=4)
          print,'parabolic model'
          LB = l / sqrt(mr-1.0)                         
          b = b0 * (1.0 + si^n / LB^n)
          dbds = n * si^(n-1.0) / ((LB^n) + (si^n)) 
       END

       2: BEGIN  
          ;; Exponential model: 
          ;; B(s) = B0 * exp(log(mr)*(s/l)^n)
          ;; where:
          ;; L: half loop length (cm)
          ;; B0: magnetic field at looptop B(L/2)=B0
          ;; mr: mirror ratio B(s)/B0
          ;; n: higher values cause stronger convergence near the footpoints
          print,'exponential model'
          dbds = alog(mr) / (l^n) * n * si^(n-1.0) ; // d (ln B) / ds
          b = b0 * exp(alog(mr) * (si/l)^n)        ;           // b(s)
       END 
    ENDCASE
    
    b[0]=b[1]
    b[n_elements(b)-1]=b[n_elements(b)-2]
    dbds[0]=0.0d
    dbds[n_elements(dbds)-1]=0.0d

    return, b

 ENDELSE 

END 

;; fp,f,e,u,s,t,b,np,ned,nud,nsd,ntd
;; rddevbin,f,e,u,s,t,b,np,ned,nud,nsd,ntd
 
;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRO rddev,f,e,u,s,t,b,np,ned,nud,nsd,ntd,file=file

;; IDL procedure to read fp.out
time0=systime(/sec)
IF ~keyword_set(file) THEN file='fp.out'
openr,unit,file,/get_lun

ned=0
nud=0
nsd=0
ntd=0

readf,unit,ned,nud,nsd,ntd ;; READ GRID SIZES

e = dblarr(ned)
u = dblarr(nud)
s = dblarr(nsd)
b = dblarr(nsd)
dbds = dblarr(nsd)
np = dblarr(nsd)

readf,unit,e ;; READ ENERGY ARRA
readf,unit,u ;; READ PITCHANGLE ARRA
readf,unit,s ;; READ POSITION ARRAY
readf,unit,b ;; READ b(s) ARRAY
readf,unit,dbds ;; READ dbds(s) ARRAY
readf,unit,np ;; READ np(s) ARRAY

readf,unit,delog,du,ds,dt ;; READ GRID RESOLUTION

i=0
n=0
tt=0

f = dblarr(ned,nud,nsd,ntd)
tmp = dblarr(nsd,nud,ned)

WHILE ~eof(unit) DO BEGIN 
   readf,unit,ti
   tt=[tt,ti]
   readf,unit,tmp
   f[0,0,0,i]=transpose(tmp)
   i++
ENDWHILE 

free_lun,unit
t=tt[1:*]
nttd=n_elements(tt)
print,'reading time ',systime(/sec)-time0

END 
;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PRO devplot,res

 f=res.f
 e=res.e
 u=res.u
 s=res.s
 t=res.t

sizes=size(f,/DIM)
ned=sizes[0]
nud=sizes[1]
nsd=sizes[2]
ntd=sizes[3]

col=interpol([10,240],ntd+1)
device,dec=0
loadct,39,/silent
window,0,xsiz=700,ysiz=690
pcs = !p.charsize
!p.charsize = 2

!p.multi=[0,2,3]
ytit='f(E,u,s,t)'

j=0
pos=nsd/2
tit='E ='+string(e[j],format='(i5)')+'keV '+'| s ='+string(s[pos]/1e9,format='(f6.1)')+'10!U9!ncm'
xtit='u = cos('+greek('phi')+')'
yr=[min(f[j,*,pos,*],/nan),max(f[j,*,pos,*],/nan)]
plot,u,f[j,*,pos,0],/nodata,/yl,tit=tit,xtit=xtit,ytit=ytit,yr=yr
FOR i=0,ntd-1 DO BEGIN 
   oplot,u,f[j,*,pos,i],ps=-4,col=col[i],lines=1
ENDFOR 

mu=nud/2
tit='E ='+string(e[j],format='(i5)')+'keV '+'| '+greek('phi')+' = '+string(acos(u[mu])/!dtor,format='(f6.1)')
xtit='position s [x10!u9!ncm]'
yr=[min(f[j,mu,*,*],/nan),max(f[j,mu,*,*],/nan)]
plot,s/1e9,f[j,mu,*,0],/nodata,/yl,tit=tit,xtit=xtit,ytit=ytit,yr=yr
FOR i=0,ntd-1 DO BEGIN 
   oplot,s/1e9,f[j,mu,*,i],ps=-2,col=col[i],lines=1
ENDFOR 

j=ned-1
tit='E='+string(e[j],format='(i5)')+'keV '+' | s ='+string(s[pos]/1e9,format='(f6.1)')+'x10!U9!ncm'
xtit='u = cos('+greek('phi')+')'
yr=[min(f[j,*,pos,*],/nan),max(f[j,*,pos,*],/nan)]
plot,u,f[j,*,pos,0],/nodata,/yl,tit=tit,xtit=xtit,ytit=ytit,yr=yr
FOR i=0,ntd-1 DO BEGIN 
   oplot,u,f[j,*,pos,i],ps=-4,col=col[i],lines=1
ENDFOR 

tit='E='+string(e[j],format='(i5)')+'keV '+' | '+greek('phi')+' = '+string(acos(u[mu])/!dtor,format='(f6.1)')+'!Uo!n'
xtit='position s [x10!u9!ncm]'
yr=[min(f[j,mu,*,*],/nan),max(f[j,mu,*,*],/nan)]
plot,s/1e9,f[j,mu,*,0],/nodata,/yl,tit=tit,xtit=xtit,ytit=ytit,yr=yr
FOR i=0,ntd-1 DO BEGIN 
   oplot,s/1e9,f[j,mu,*,i],ps=-2,col=col[i],lines=1
ENDFOR

tit=greek('phi')+' = '+string(acos(u[mu])/!dtor,format='(f6.1)')+'!Uo!n'+' | s = '+string(s[pos]/1e9,format='(f6.1)')+'x10!U9!ncm'
xtit='energy [keV]'
yr=[min(f[*,mu,pos,*],/nan),max(f[*,mu,pos,*],/nan)]
plot,e,f[*,mu,pos,0],/xl,/yl,/nodata,tit=tit,xtit=xtit,ytit=ytit,yr=yr
FOR i=0,ntd-1 DO BEGIN 
   oplot,e,f[*,mu,pos,i],ps=0,col=col[i],lines=0
ENDFOR

pos=nsd-2
mu=nud-2
tit=greek('phi')+' = '+string(acos(u[mu])/!dtor,format='(f6.1)')+'!Uo!n'+' | s = '+string(s[pos]/1e9,format='(f6.1)')+'x10!U9!ncm'
xtit='energy [keV]'
plot,e,f[*,mu,pos,0],/xl,/yl,/nodata,tit=tit,xtit=xtit,ytit=ytit,yr=[min(f[*,mu,pos,*],/nan),max(f[*,mu,pos,*],/nan)]
FOR i=0,ntd-1 DO BEGIN 
   oplot,e,f[*,mu,pos,i],ps=0,col=col[i],lines=0
ENDFOR

!p.charsize=pcs
device,dec=1
!p.multi=0
END 

;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRO rddevbin,result,file=file

 IF ~keyword_set(file) THEN file='fp.out'
 openr,unit,file,/get_lun
 
 grid=lonarr(4)
 res=dblarr(4)

 readu,unit,grid
 readu,unit,res

 narray=grid[0]*grid[1]*grid[2]
 e=dblarr(grid[0])
 u=dblarr(grid[1])
 s=dblarr(grid[2])
 b=dblarr(grid[2])
 dbds=dblarr(grid[2])
 np=dblarr(grid[2]) 
 nt=0l

 readu,unit,e
 e*=510.99891
 readu,unit,u
 readu,unit,s
 readu,unit,b
 readu,unit,dbds
 readu,unit,np

 readu,unit,nt

 t=dblarr(nt)
 d=dblarr(nt*narray)

 readu,unit,t
 readu,unit,d
 free_lun,unit

 grid[3]=nt
 f=dblarr(grid)
 i=lindgen(narray)
 FOR j=0,nt-1 DO $
  f[0:grid[0]-1,0:grid[1]-1,0:grid[2]-1,j] $
  =transpose(reform(d[i+narray*j],grid[2],grid[1],grid[0]))

 result={f:f,e:e,u:u,s:s,t:t,b:b,np:np,ned:grid[0] $
         ,nud:grid[1],nsd:grid[2],ntd:grid[3],dbds:dbds}

END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO writebfield,s,dbds

 openw,unit,'bfielddata.txt',/get_lun
 printf,unit,n_elements(dbds)
 printf,unit,s
 printf,unit,dbds
 free_lun,unit

END 


PRO fp,p,res,bfield=bfield,plot=plot,file=file,dontrun=dontrun

 IF n_elements(p) EQ 0 THEN BEGIN  

    ;; grid sizes
    ned=30
    nud=80
    nsd=80
    ntd=100
    ;; energy limits
    emin=100.0
    emax=200.0
    ;; position limits
    smin=-2d9
    smax=2d9
    ;; magnetic field model
    mr=2d
    nexp=2d
    bmodel=0
    ;; plasma density model
    nmin=1d10
    nmax=1d10
    ;; max time of simulation
    tmax=5.0
    ;; injection function
    d0=3.0
    dp=2.3
    p0=1.0
    dx=1e8
    x0=0.0
    t0=2.0 ;;0.25
    tau=0.0 ;;0.1
    lambda=1d20

    ;; to set a continuous constant injection: tau=0, t0 ne 0, run[4]=1
    run0=1
    run1=0
    run2=1
    run3=1
    run4=0
    run5=0

 ENDIF ELSE BEGIN 

    ;; grid sizes
    ned=p[0]
    nud=p[1]
    nsd=p[2]
    ntd=p[3]
    ;; energy limits
    emin=p[4]
    emax=p[5]
    ;; position limits
    smin=p[6]
    smax=p[7]
    ;; magnetic field model
    mr=p[8]
    nexp=p[9]
    bmodel=p[10]
    ;; plasma density model
    nmin=p[11]
    nmax=p[12]
    ;; max time of simulation
    tmax=p[13]
    ;; injection function
    d0=p[14]
    dp=p[15]
    p0=p[16]
    dx=p[17]
    x0=p[18]
    t0=p[19] ;;0.25
    tau=p[20] ;;0.1
    lambda=p[21]

    ;; to set a continuous constant injection: tau=0, t0 ne 0, run[4]=1
    run0=p[22]
    run1=p[23]
    run2=p[24]
    run3=p[25]
    run4=p[26]
    run5=p[27]

 ENDELSE 

 run=[run0,run1,run2,run3,run4,run5]
 arguments=''
 ;; dt override
;dt=1e-4
;arguments=' -t '+strtrim(string(dt),2)

 ;; prepare to write
 grid=fix([ned,nud,nsd,ntd])
 par1=[emin,emax,smin,smax,mr,nexp]
 par2=fix(bmodel)
 par3=[nmin,nmax,tmax,d0,dp,p0,dx,x0,t0,tau*run[5],lambda]

 ;; write input file
 openw,unit,'fp.in',/get_lun
 printf,unit, fix(grid)
 printf,unit, par1
 printf,unit, par2
 printf,unit, par3
 printf,unit, fix(run)
 free_lun,unit

 IF (bmodel EQ 2) THEN BEGIN 
    ;; restore,'bfielddata.sav' ;bfield[100,2]
    if keyword_set(bfield) then BEGIN 
       s_field=bfield[*,0]
       dbds_field=deriv(s_field,alog(bfield[*,1]))
    ENDIF ELSE BEGIN
       s_field=interpol([smin,smax],nsd) 
       b=bfield(s_field,mr=mr,n=nexp,/prec,dbds=dbds_field)
       bfield=[[s_field],[b]]
    ENDELSE 
    writebfield,s_field,dbds_field
    ;;plot,s_field,dbds_field
 ENDIF 

 time0=systime(/sec)
 command='~/solar/fp/dev'

if ~keyword_set(dontrun) then BEGIN 
   print,'Running Fokker-Planck.'
   spawn, command+arguments
   ctime=systime(/sec)-time0

   ;; read output
   ;;rddev,f,e,u,s,t,b,np,ned,nud,nsd,ntd
   rddevbin,res
   
   IF (bmodel EQ 2) THEN res.b=interpol(bfield[*,1],res.nsd,/spline)
   
   ;; plotting
   if keyword_set(plot) then devplot,res
ENDIF ELSE print,'fp.in written'

END
