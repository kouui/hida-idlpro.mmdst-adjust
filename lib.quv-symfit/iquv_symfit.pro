; LIBRARY FOR FITTING SYMMETRY Q/U/V PROFILE OF SUNSPOT
;
; HISTORY : 
; 2021.02.13  U.K.  CREATED, WORKS

;; EXTERNAL DEPENDENCY
@mmdst_lib

;------------------------------------------------------------------------
;; FIND THE POSITION IN AN ARRAY 
;; WITH A VALUE MOSTLY CLOSE TO THE TARGET VALUE
FUNCTION ARRAY1D_FIND_NEAREST, array, target, value=value

    ss = size(array)
    if ss[0] ne 1 then throw_error, "ndim of array !=1"
    
    value = min( abs(array-target), pos )
    
    return, pos
END


;------------------------------------------------------------------------
;; CALCULATE THE RANGE OF OBSORPTION LINE AND CONTINUUM
;; - FOR OBSORPTION LINE : FIT THE MINIMUM AND INTERPOLATE THE GRID
;; - FOR CONTINUUM : INTERPOLATE TO HAVE THE SAME NUMBER OF POINTS WITH LINE
;; xl1, xl2, xc1, xc2 : in pixel unit
FUNCTION CALCULATE_FIT_RANGE, xl1, xl2, xc1, xc2, sp, xl3=xl3, xc3=xc3, xls=xls, xcs=xcs

    xl1   = fix(xl1)
    xl2   = fix(xl2)
    xc1   = fix(xc1)
    xc2   = fix(xc2)

    xminv = min(sp[xl1:xl2], xmin)
    ;; fit minimum to find center
    x0 = x0fit(sp,xmin+xl1,8) & x0 = x0[0]
    nx = MINVAL(fix( abs(x0-xl1) ), fix( abs(x0-xl2) ))
    nxl = nx
    npoints = 2*nx+1
    xl3 = [x0-nx, x0, x0+nx]

    xc0 = (xc1 + xc2) / 2
    nx = fix( abs(xc0-xc1) )
    xc3 = [xc0-nx, xc0, xc0+nx]

    xls = indgen(npoints, start=xl3[0], /float)
    ;; shrink the interval in continuum grid
    ;; to match 
    ;; - the given range
    ;; - npoint in line grid
    xcs = indgen(npoints, start=xc3[0], /float, increment=float(nx)/float(nxl))

    return, npoints
END

;------------------------------------------------------------------------
;; GIVEN THE CALCULATED GRID
;; INTERPOLATE IQUV PROFILE
;; s4, (nx,ny,4,1), 4 : I,Q,U,V, 1 : nseq
;; spil, (nx,npoints,4,1), 4 : I,Q,U,V, 1 : nseq
;; spic, (nx,npoints,4,1), 4 : I,Q,U,V, 1 : nseq
FUNCTION INTERP_PROFILE, s4, xls, xcs, spil=spil, spic=spic

    ss = size(s4)
    if  ss[0] ne 4 then throw_error, "ndim of s4 !=4"
    if  ss[3] ne 4 then throw_error, "size of the 3rd dimension of s4 !=4"
    if  ss[4] ne 1 then throw_error, "size of the 4th dimension of s4 !=1"
    nx = ss[1]
    ny = ss[2]
    xdata = indgen(ny)

    ;; interpolate line
    npoints = n_elements(xls)
    spil = fltarr(nx,npoints,4,1)
    for i=0,3 do begin
    for j=0,nx-1 do begin
        ydata = reform(s4[j,*,i,0])
        spil[j,*,i,0] = interpolate(ydata,xls,/cubic)
    endfor
    endfor
    ;; interpolate continuum
    npoints = n_elements(xcs)
    spic = fltarr(nx,npoints,4,1)
    for i=0,3 do begin
    for j=0,nx-1 do begin
        ydata = reform(s4[j,*,i,0])
        spic[j,*,i,0] = interpolate(ydata,xcs) ;; bilinear
    endfor
    endfor

    return, 1b

END

;------------------------------------------------------------------------
;; spil, (nx,npoints,4,nseq), 4 : I,Q,U,V
;; arr, (nx,npoints,3,2,nseq), 3 : Q,U,V, 2 : line,continuum
FUNCTION MAKE_HALF_ARRAY, spil
    
    ss = size(spil)
    if ss[0] ne 4 then throw_error, "ndim of spil !=4"
    npoints = ss[2]
    npoints_half = (npoints-1) / 2
    
    arr = fltarr(ss[1],npoints_half,3,2,ss[4])

    return, arr
END
;------------------------------------------------------------------------
;; error array given spil
FUNCTION MAKE_ERROR_ARRAY, spil, err
    
    arr = MAKE_HALF_ARRAY(spil)
    arr[*,*,*,*,*] = 1.0 * err

    return, arr
END

;------------------------------------------------------------------------
;; half arr, (nx,npoints,3,2,nseq), 3 : Q,U,V, 2 : line,continuum
;; spil, (nx,npoints,4,nseq), 4 : I,Q,U,V
;; weight array given spil
;; nx range :
;; active region [0:xaq-1]
;; quiescent region [xaq:nx-1]
FUNCTION MAKE_WEIGHT_ARRAY, spil, xaq, $
            wtQ=wtQ, wtU=wtU, wtV=wtV, $
            wtL=wtL, wtC=wtC, wt_AR=wtAR, wt_QR=wtQR
    
    if not keyword_set(wtQ) then wtQ = 1.0
    if not keyword_set(wtU) then wtU = 1.0
    if not keyword_set(wtV) then wtV = 1.0
    if not keyword_set(wtL) then wtL = 1.0
    if not keyword_set(wtC) then wtC = 1.0
    if not keyword_set(wt_AR) then wt_AR = 1.0
    if not keyword_set(wt_QR) then wt_QR = 1.0


    arr = MAKE_HALF_ARRAY(spil)
    ss = size(arr) & nx = ss[1]
    
    arr[*,*,*,*,*] = 1.0
    arr[*,*,0,*,*] *= wtQ
    arr[*,*,1,*,*] *= wtU
    arr[*,*,2,*,*] *= wtV
    arr[*,*,*,0,*] *= wtL
    arr[*,*,*,1,*] *= wtC
    arr[0:xaq-1,*,*,*,*] *= wt_AR
    arr[xaq:nx-1,*,*,*,*] *= wt_QR

    return, arr
END

;------------------------------------------------------------------------
;; result/target array given spil and spic
;; half arr, (nx,npoints,3,2,nseq), 3 : Q,U,V, 2 : line,continuum
;; spil, (nx,npoints,4,nseq), 4 : I,Q,U,V
;; nx range :
;; active region [0:xaq-1]
;; quiescent region [xaq:nx-1]
FUNCTION MAKE_RESULT_ARRAY, spil, spic, xaq, target=target, zero_netV = zero_netV, $
                neglect_QRV=neglect_QRV, neglect_continuum=neglect_continuum

    arr = MAKE_HALF_ARRAY(spil)
    ;; diff between left and right of 
    ;; a (anti-)symmetry profile should be zero
    if keyword_set(target) then return, arr

    ss = size(arr) & nph = ss[2] & nx = ss[1] & nseq = ss[5]
    ;; Line AR : Q, U, symmetry
    for i=0,1 do begin
        arr[0:xaq-1,*,i,0,*] = spil[0:xaq-1,0:nph-1,i+1,*] - reverse(spil[0:xaq-1,nph+1:2*nph,i+1,*],2)
    endfor
    ;; Line AR : V, anti-symmetry
    arr[0:xaq-1,*,2,0,*] = spil[0:xaq-1,0:nph-1,3,*] + reverse(spil[0:xaq-1,nph+1:2*nph,3,*],2)
    if not keyword_set(zero_netV) then arr[0:xaq-1,*,2,0,*] -= spil[0:xaq-1,nph,3,*]

    ;; Line QR : Q, U, V zero/flat
    for i=0,2 do begin
        arr[xaq:nx-1,*,i,0,*] = 0.5*(spil[xaq:nx-1,0:nph-1,i+1,*] + reverse(spil[xaq:nx-1,nph+1:2*nph,i+1,*],2))
        if not keyword_set(zero_netV) then arr[xaq:nx-1,*,i,0,*] -= rebin(arr[xaq:nx-1,*,i,0,*],nx-xaq,1,1,1,nseq)
    endfor
    if keyword_set(neglect_QRV) then arr[xaq:nx-1,*,2,0,*] = 0


    ;; Continuum : zero/flat
    for i=0,2 do begin
        arr[*,*,i,1,*] = 0.5*(spic[*,0:nph-1,i+1,*] + spic[*,nph+1:2*nph,i+1,*])
        if not keyword_set(zero_netV) then arr[*,*,i,1,*] -= rebin(arr[*,*,i,1,*],1,nph,1,1,nseq)
    endfor
    if keyword_set(neglect_continuum) then arr[*,*,0:1,1,*] = 0

    return, arr
END

;------------------------------------------------------------------------
;; parinfo struct for mpfit
; create parinfo for a single parameter in mpfit
FUNCTION INIT_PARINFO 
    parinfo = {parinfo, value:0.0, fixed:0b, limited:[1b,1b], limits:[0.0,1.0], parname:'sc'}
    return, parinfo
END

; create parinfo for 5 parameters in mpfit
FUNCTION MAKE_PARINFO, xn, tn, xc, tc, sc, th_vs, conf_symfit
    npar = 6
    parinfo = [] & for i=0,npar-1 do parinfo = [parinfo,init_parinfo()]
    ;; xn
    i=0
    parinfo[i].value = xn
    parinfo[i].parname = 'xn'
    parinfo[i].limits = [-0.1,+0.1]
    ;; tn
    i=1
    parinfo[i].value = tn
    parinfo[i].parname = 'tn'
    parinfo[i].limits = [-90.0*!dtor,+90.0*!dtor]
    ;; xc
    i=2
    parinfo[i].value = xc
    parinfo[i].parname = 'xc'
    parinfo[i].limits = [-0.1,+0.1]
    ;; tc
    i=3
    parinfo[i].value = tc
    parinfo[i].parname = 'tc'
    parinfo[i].limits = [-90.0*!dtor,+90.0*!dtor]
    ;; sc
    i=4
    parinfo[i].value = sc
    parinfo[i].parname = 'sc'
    parinfo[i].limits = [-0.3,+0.3]
    parinfo[i].fixed = conf_symfit.fix_sc
    ;; th_vs
    i=5
    parinfo[i].value = th_vs
    parinfo[i].parname = 'th_vs'
    parinfo[i].limits = [-45.0,+45.0]
    parinfo[i].fixed = conf_symfit.fix_th_vs
    return,parinfo
END

;------------------------------------------------------------------------
;; spils, (nx,npoints,4,nseq), 4 : I,Q,U,V
;; dst3c, (nseq), array of struct
FUNCTION MODEL_PARDST, cx, par
    common symfit, spils, spics, dst3c, xaq

    xn=par[0] & tn=par[1] & xc=par[2] & tc=par[3] & sc=par[4] & th_vs=par[5]

    ss = size(spils) & nx = ss[1] & np=ss[2] & nst=ss[3] & ns=ss[4]
    nph = (np-1) / 2
    spmls = fltarr(nx,np,4,ns)
    spmcs = fltarr(nx,np,4,ns)
    if ns eq 1 then begin
        spmls = reform(spmls,nx,np,4,ns,/overwrite)
        spmcs = reform(spmcs,nx,np,4,ns,/overwrite)
    endif

    for i=0,ns-1 do begin
        dst = dst3c[i]
        mm = UPDATE_MMDST(dst, xn, tn, xc, tc, sc, th_vs=th_vs)
        spmls[*,*,*,i] = UPDATE_S3(reform(spils[*,*,*,i],nx,np,nst), mm)
        spmcs[*,*,*,i] = UPDATE_S3(reform(spics[*,*,*,i],nx,np,nst), mm)
    endfor
    res = MAKE_RESULT_ARRAY(spmls, spmcs, xaq, /zero_netV, /neglect_QRV, /neglect_continuum)


    print, "[SymFit] mean of Absolute error : ", mean(res)
    return, res
END

;------------------------------------------------------------------------
;; s4, (nx, ny, 4, nseq)
;; dst3, (3, 11, nseq)
FUNCTION MPFIT_PARDST, conf_symfit, s4, dst3, parinit, parinfo, err, weight, niter=niter, quiet=quiet
    common symfit, spils, spics, dst3c, xaq
    dst3c = dst3
    xaq = conf_symfit.xaq

    ss = size(s4)
    if ss[0] ne 4 then begin
    if ss[0] eq 3 then s4 = reform(s4,ss[1],ss[2],ss[3],1,/overwrite) else throw_error, "ndim of s4 == ",ss[0]
    endif
    ss = size(s4)
    ns = ss[4] & nx = ss[1] & ny = ss[2] & nst = ss[3]
    sp = reform(rebin(s4[xaq:nx-1,*,0,0],1,ny,1,1))
    ;x0 = nx/2
    ;ret = CALCULATE_FIT_RANGE(460, 445, 120, 110, reform(s4[0,0,*,0]),xls=xls, xcs=xcs)
    ret = CALCULATE_FIT_RANGE(conf_symfit.ixL1, conf_symfit.ixL2, conf_symfit.ixC1, conf_symfit.ixC2,sp,xls=xls, xcs=xcs)
    ;ret = CALCULATE_FIT_RANGE(460, 445, 277, 267, reform(s4[0,0,*,0]),xls=xls, xcs=xcs)
    np = n_elements(xls)
    ;nph = (np-1)/2
    spils = fltarr(nx,np,4,ns)
    spics = fltarr(nx,np,4,ns)
    if ns eq 1 then begin
        spils = reform(spils,nx,np,4,ns,/overwrite)
        spics = reform(spics,nx,np,4,ns,/overwrite)
    endif
    for i=0,ns-1 do begin
        ret = INTERP_PROFILE(reform(s4[*,*,*,i],nx,ny,nst,1), xls, xcs, spil=spil, spic=spic)
        spils[*,*,*,i] = spil[*,*,*,*]
        spics[*,*,*,i] = spic[*,*,*,*]
    endfor
    ;; display and check interpolated obsorption line
    showiquvr,reform(spils[*,*,*,0],nx,np,4), bin=1, wid=17
    ;stop
    

    weights = MAKE_WEIGHT_ARRAY(spils,xaq,wtQ=weight.wtQ,wtU=weight.wtU,wtV=weight.wtV,$
                wtL=weight.wtL,wtC=weight.wtC,wt_AR=weight.wtAR,wt_QR=weight.wtQR)
    
    if not keyword_set(niter) then niter=20
    if not keyword_set(quiet) then quiet=1b

    errs =  MAKE_ERROR_ARRAY(spils, err)
    targets = MAKE_RESULT_ARRAY(spils, spics, xaq, /target)
    ;;stop
    ret = mpfitfun('MODEL_PARDST',0,targets,errs,parinit,weight=weights,errmsg=errmsg, $
				quiet=quiet,maxiter=niter,niter=k,yfit=sfit,bestnorm=eval1,parinfo=parinfo)


    return,ret

END


sconf = {symfit_struct,$
    ixL1 : 430, $
    ixL2 : 460, $
    ixC1 : 100, $
    ixC2 : 130, $
    iyA1 : 70, $
    iyA2 : 150, $
    iyQ1 : 180, $
    iyQ2 : 240, $
    xaq  : 80, $
    error: 2E-4, $
    niter: 60, $
    fix_sc : 1b, $
    fix_th_vs : 1b, $
    corr_Icrtk: 1b $
    }
weight = {weight_struct, wtQ:1.0, wtU:1.0, wtV:1.0, wtL:1.0, wtC:1.0, wtAR: 1.0, wtQR:1.0}
parinit = [0.00300, -0.06286, 0.02500, 0.12566, 0.0, 0.0]
parinfo = MAKE_PARINFO(parinit[0], parinit[1], parinit[2], parinit[3], parinit[4], parinit[5], sconf)

file = '/nwork/kouui/data-lvl1/dstpol/20220110.giant-prominence/spec/cal/sav.for-mmdst-adjust'
file = file + '/step5.s.ar.sav'
restore, file
dst3 = [dst]
ss = size(s) & nx = ss[1] & ny = ss[2] & nst = ss[3]-1


s4d = fltarr(sconf.iyA2-sconf.iyA1+sconf.iyQ2-sconf.iyQ1, ny, nst)
ss = size(s4d) & nx = ss[1] & ny = ss[2] & nst = ss[3]
s4d = reform(s4d, nx, ny, nst, 1)
s4d[0:sconf.xaq-1,*,*,0] = s[sconf.iyA1:sconf.iyA2-1,*,0:3]
s4d[sconf.xaq:nx-1,*,*,0] = s[sconf.iyQ1:sconf.iyQ2-1,*,0:3]

;showiquvr, s[*,*,0:3]
;draw, (240+(2+nx)*2)*[1,1], [2,2+ny], color=150, thick=2, line=2
showiquvr, reform(s4d), wid=2
draw, [0,3*2+nx*4], sconf.ixL1*[1,1], color=150, thick=2, line=2
draw, [0,3*2+nx*4], sconf.ixL2*[1,1], color=150, thick=2, line=2
draw, [0,3*2+nx*4], sconf.ixC1*[1,1], color=150, thick=2, line=2
draw, [0,3*2+nx*4], sconf.ixC2*[1,1], color=150, thick=2, line=2

pars_fit = MPFIT_PARDST(sconf, s4d, dst3, parinit, parinfo, sconf.error, weight, niter=sconf.niter)

mm = update_mmdst(dst, pars_fit[0], pars_fit[1], pars_fit[2], pars_fit[3], pars_fit[4],th_vs=pars_fit[5])
sr = UPDATE_S3(s[*,*,0:3],mm)
showiquvr, sr, wid=3

print, pars_fit

END