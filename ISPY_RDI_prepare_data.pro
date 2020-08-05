pro ISPY_RDI_prepare_data

lambda = 3.8d-6 ;needed to compute the distance of the waffle spots
plsc = 27.19  ;mas/px

datapath = '/data/beegfs/astro-storage/groups/launhardt/amueller/NACO/DACE/data_raw/'
; datapath = '/home/amueller/Downloads/ISPY_RDI/data/'

resdir = datapath+'../Lp/'
spawn, 'mkdir -p '+resdir

cut = 30. ;cutting radius

;images are 50x50px

; dim = 50.
; xa = dindgen(dim*dim)
; errdum = dblarr(dim^2.)
; errdum[*] = 1.

; see line 45
; mask_t = shift(dist(dim), dim/2., dim/2.)
; mask = mask_t ge 3. and mask_t le 15.
; idxmask = where(mask ne 0.)
; npxmask = double(n_elements(idxmask))

;-------------------------------------------------------------------------

file1 = file_search(datapath+'*master.fits', count=n)
file2 = file_search(datapath+'*PA.rdb', count=n2)

for i=0,n-1 do begin
;for i=107,107 do begin

  print, i, ' ', n, ' ', file1[i]
  sz = (file_info(file1[i])).size

  if (sz gt 0.1) then begin

		readcol, file2[i], paral, format='d', /silent
;     paral = mrdfits(file2[i],0,/silent)
;     paral = -1.*paral

    sz = 0
    hdr = headfits(file1[i], exten=0, /silent)
    date = get_eso_keyword(hdr, 'DATE-OBS')
    object = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO OBS TARG NAME'),/rem)
    naxis3 = get_eso_keyword(hdr, 'NAXIS3')
    object_orig = object
    date = strmid(date,0,19)
    date = repstr(date, ':', '_')
    object = strtrim(object,2)
    object = repstr(object, ' ', '_')
    object = repstr(object, '__', '_')
    object = repstr(object, '___', '_')
    object = repstr(object, '____', '_')
    object = repstr(object, '_____', '_')
    print, object_orig, ' ', object, ' ', date

;     if (strmatch(object, '*HIP0*') eq 1 or strmatch(object, '*HIP1*') eq 1 or strmatch(object, '*HIP2*') eq 1 or strmatch(object, '*HIP3*') eq 1 or strmatch(object, '*HIP4*') eq 1 or strmatch(object, '*HIP5*') eq 1 or strmatch(object, '*HIP6*') eq 1 or strmatch(object, '*HIP7*') eq 1 or strmatch(object, '*HIP8*') eq 1 or strmatch(object, '*HIP9*') eq 1) then begin
; 
;       pos = strpos(object, 'HIP')
;       num = strmid(object, 3, strlen(object)-3)
;       object = 'HIP_'+num
; 
;     endif
;     
; 
;     if (strmatch(object, '*HD0*') eq 1 or strmatch(object, '*HD1*') eq 1 or strmatch(object, '*HD2*') eq 1 or strmatch(object, '*HD3*') eq 1 or strmatch(object, '*HD4*') eq 1 or strmatch(object, '*HD5*') eq 1 or strmatch(object, '*HD6*') eq 1 or strmatch(object, '*HD7*') eq 1 or strmatch(object, '*HD8*') eq 1 or strmatch(object, '*HD9*') eq 1) then begin
; 
;       pos = strpos(object, 'HD')
;       num = strmid(object, 2, strlen(object)-2)
;       object = 'HD_'+num
; 
;     endif
; 
;     if (strmatch(object, '*TYC0*') eq 1 or strmatch(object, '*TYC1*') eq 1 or strmatch(object, '*TYC2*') eq 1 or strmatch(object, '*TYC3*') eq 1 or strmatch(object, '*TYC4*') eq 1 or strmatch(object, '*TYC5*') eq 1 or strmatch(object, '*TYC6*') eq 1 or strmatch(object, '*TYC7*') eq 1 or strmatch(object, '*TYC8*') eq 1 or strmatch(object, '*TYC9*') eq 1) then begin
; 
;       pos = strpos(object, 'TYC')
;       num = strmid(object, 3, strlen(object)-3)
;       object = 'TYC_'+num
; 
;     endif
; 
;     if (strmatch(object, '*No_name*') eq 1) then begin
; 
;       object = ''
;       print, file1[i]
;       read, 'Provide object name: ', object
; 
;     endif


    d = mrdfits(file1[i], 0, hdr, /silent)
    ;idxo = where(strmatch(hdr, 'OBJECT*'))
    ;hdr[idxo] = "OBJECT  = '"+object+"'      / Original target."

    ;check if number of frames is equal to the number of angles
    if (n_elements(d[0,0,*]) ne n_elements(paral)) then begin

      print, 'Number of frames different from number of angles! Stop.'
      stop

    endif

    dim = n_elements(d[*,0,0])
    d1 = d[dim/2-cut:dim/2+cut-1,dim/2-cut:dim/2+cut-1,*,0]
    ;d2 = d[dim/2-50:dim/2+49,dim/2-50:dim/2+49,*,1]

    nframes = n_elements(d1[0,0,*])

     ;check if a frame has NANs
    flag = intarr(nframes)
    for j=0,nframes-1 do flag[j] = total(d1[*,*,j])
    idx = where(flag ne 0)

    d1 = d1[*,*,idx]
    paral = paral[idx]
    nframes = n_elements(idx)
    naxis3 = nframes


;     medim1 = median(d1, dim=3)
;     sim1 = dblarr(nframes)
;     for j=0,nframes-1 do sim1[j] = stddev(d1[*,*,j]-medim1)
;     resistant_mean, sim1, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad
; 
;     d1 = d1[*,*,idxgood]
;     d2 = d2[*,*,idxgood]
    
;     paral = paral[idxgood]

    ;normalize data
;     naxis3 = n_elements(idxgood)
    normfactor = dblarr(nframes)
;     normfactorH3 = dblarr(n_elements(idxgood))
    for j=0,naxis3-1 do begin
    
      normfactor[j] = total(d1[*,*,j])
;       normfactorH3[j] = total(d2[*,*,j])
      d1[*,*,j] = d1[*,*,j]/normfactor[j]
;       d2[*,*,j] = d2[*,*,j]/normfactorH3[j]
        
    endfor
    
    writefits, resdir+object+'_'+date+'_Lp.fits', d1, hdr
    writefits, resdir+object+'_'+date+'_Lp.fits', paral, /append
;     writefits, resdir+object+'_'+date+'_H3.fits', d2, hdr
;     writefits, resdir+object+'_'+date+'_H3.fits', paral, /append

  ;  spawn, 'mkdir -p '+resdir+object+'_'+date
  ;  writefits, resdir+object+'_'+date+'/img_H2_dc.fits', d1, hdr
  ;  writefits, resdir+object+'_'+date+'/img_H3_dc.fits', d2, hdr
  ;  spawn, 'cp '+file2[i]+' '+resdir+object+'_'+date+'/vec_H2_paral.fits'
  ;  spawn, 'cp '+file2[i]+' '+resdir+object+'_'+date+'/vec_H3_paral.fits'


    dim = n_elements(d1[*,0,0])
    flux = dblarr(naxis3)
    sdevflux = dblarr(naxis3)
    meanflux = dblarr(naxis3)
;     fluxH3 = dblarr(naxis3)
;     sdevfluxH3 = dblarr(naxis3)
;     meanfluxH3 = dblarr(naxis3)
  
    mask_t = shift(dist(dim), dim/2., dim/2.)
    mask1 = mask_t ge 3.  ;avoid the very center
    dspot = 25.;(14.*lambda/8./!dtor)*3600.*1000./plsc - 10. 
    mask2 = mask_t le dspot ;avoid potential waffle pattern
    mask = mask1*mask2

    idxmask = where(mask ne 0.)
    npxmask = double(n_elements(idxmask))

    im1d = dblarr(naxis3, npxmask)
;     im1dH3 = dblarr(naxis3, npxmask)
    
    for j=0,naxis3-1 do begin
        
      t1 = (d1[*,*,j])[*]
      im1d[j,*] = t1[idxmask]
      flux[j] = total(im1d[j,*])
      sdevflux[j] = stddev(im1d[j,*])
      meanflux[j] = mean(im1d[j,*])
      
;       t1 = (d2[*,*,j])[*]
;       im1dH3[j,*] = t1[idxmask]
;       fluxH3[j] = total(im1dH3[j,*])
;       sdevfluxH3[j] = stddev(im1dH3[j,*])
;       meanfluxH3[j] = mean(im1dH3[j,*])
      
    endfor
    
    target = object
    
;     im1d = im1dH2
;     flux = fluxH2
;     sdevflux = sdevfluxH2
;     meanflux = meanfluxH2
;     normfactor = normfactorH2
    fn = resdir+object+'_'+date+'_Lp.sav'
;     fn_all[i] = fn
    save, target, im1d, naxis3, paral, flux, sdevflux, meanflux, mask, npxmask, idxmask, normfactor, filename=fn
    
;     im1d = im1dH3
;     flux = fluxH3
;     sdevflux = sdevfluxH3
;     meanflux = meanfluxH3
;     normfactor = normfactorH3
;     fn = resdir+object+'_'+date+'_H3.sav'
;     save, target, im1d, naxis3, paral, flux, sdevflux, meanflux, mask, npxmask, idxmask, normfactor, filename=fn    
;     

  endif else begin

    count0 = count0+1

  endelse


endfor


;-------------------------------------------------------------------------

stop
end


