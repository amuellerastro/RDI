@dist.pro
@readcol.pro
@remchar.pro
@gettok.pro
@strnumber.pro
@proceeding_text.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@fxmove.pro
@get_eso_keyword.pro
@rot.pro
@writefits.pro
@mkhdr.pro
@sxaddpar.pro
@get_date.pro
@daycnv.pro
@sxdelpar.pro
@sxpar.pro
@check_fits.pro
@fxaddpar.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@host_to_ieee.pro
@fits_ascii_encode.pro
@resistant_mean.pro
@poly.pro
@strsplit.pro
@mpfitfun.pro
@mpfit.pro
@ms_ssim.pro
@curvefit.pro
@gauss_smooth.pro
@readfits.pro
@filter_image.pro
@array_indices.pro
@factor.pro
@prime.pro
@psf_gaussian.pro
@gaussian.pro
@convolve.pro
@real_part.pro

;Wang 2004 Image Quality Assessment: From Error Visibility to Structural Similarity
;http://www.cns.nyu.edu/pub/lcv/wang03-preprint.pdf

; function func_fit_psf_amplitude, x, p, _extra=st
; 
; 	t1 = st.t1
; 	t2 = st.t2
; 	
; 	fit = t2*p[0]
; 	
; 	return, fit
; 	
; end

; pro func_fit_psf_amplitude, x, A, F
; 
; common data, st
; 
; 	t1 = st.t1
; 	t2 = st.t2
; 	
; 	F = t2*A[0]
;     
; end

pro ISPY_RDI_find_ref
Compile_opt idl2

;just search for the references
;-------------------------------------------------------------------------

filter = 'Lp'
;read, 'Filter (e.g. H2): ', filter

datapath = '/data/beegfs/astro-storage/groups/launhardt/amueller/NACO/DACE/data/'+filter+'/'
reffile = '/home/amueller/work/IDLlibs/AO/NACO/ISPY_RDI_targets.dat'

nbest = 100 ;number of best frames

qmask = 'y'
read, 'Mask (y/n): ', qmask

if (qmask eq 'y') then resdir = datapath+'../Reference_'+filter+'_mask/' else resdir = datapath+'../Reference_'+filter+'_nomask/'

tmpfile = file_search(datapath+'*.fits')
hdr = headfits(tmpfile[0], exten=0, /silent)
tmpfile = ''
dim = double(get_eso_keyword(hdr, 'NAXIS3'))

if (filter eq 'Lp') then lambda = 3.8d-6
;if (filter eq 'H2') then lambda = 1.593d-6
;if (filter eq 'H3') then lambda = 1.667d-6
diam = 8.2
plsc = 27.19d-3
fwhm = lambda/diam*206265.d0/plsc

;-------------------------------------------------------------------------
spawn, 'mkdir -p '+resdir
readcol, reffile, id, systype, flag, format='a,a,i', /silent
;-------------------------------------------------------------------------

file = file_search(datapath+'*.sav', count=nfiles)
fitsfile = file_search(datapath+'*.fits', count=nfitsfiles)
ref_fitsfile = file_search(datapath+'*.fits', count=nrfitsfiles)
ref_file = file_search(datapath+'*.sav', count=nrfitsfiles)

if (nfitsfiles ne nrfitsfiles) then begin

    print, 'Number of fits files in average folder different! Stop.'
    stop

endif

if (nfiles ne nfitsfiles) then stop

;idx = (where(strmatch(fitsfile, '*HD97048_2017-05-03T02_22_34*') eq 1))[0]

;loop to go through every science target
for isci=0,nfiles-1 do begin

    t = systime(1)

    sci_im = readfits(fitsfile[isci],exten_no=0,/silent)
    file_id = strmid(fitsfile[isci], strpos(fitsfile[isci], '/', /reverse_search)+1, strlen(fitsfile[isci])-strpos(fitsfile[isci], '/', /reverse_search)-6)

    restore, file[isci]
	  sci_target = target
	  sci_im1d = im1d
    sci_naxis3 = naxis3

;
    print, 'Working on file '+strcompress(sci_target,/rem)+' '+strcompress(isci+1,/rem);+' / '+strcompress(nfiles,/rem)

	;measure following metric
	;Wang 2004 Image Quality Assessment: From Error Visibility to Structural Similarity
	;http://www.cns.nyu.edu/pub/lcv/wang03-preprint.pdf
	
; 	if keyword_set(pca) then begin
	
;     best_ssim = dblarr(sci_naxis3, nfiles)
    best_ssim = dblarr(sci_naxis3, nbest)
    diff_im = dblarr(dim, dim, sci_naxis3)
    ref_frame = dblarr(dim, dim, sci_naxis3)
    best_frames = dblarr(dim, dim, nbest)
    ssim = dblarr(sci_naxis3, nfiles, 2000) ;500 is just a dummy, usually we don't have more than 200 frames pr star
    ssim_file = strarr(sci_naxis3, nfiles, 2000)
    ssim_file_frame = intarr(sci_naxis3, nfiles, 2000)
    bestidx = intarr(sci_naxis3, nbest, 2)   ;best 100 frames
    best_ssim_file = strarr(sci_naxis3, nbest)
    best_ssim_file_frame = intarr(sci_naxis3, nbest)
	
	;loop to go through every data set to find the reference frame to be used
	for iref=0,nfiles-1 do begin

        ref_im = readfits(ref_fitsfile[iref],exten_no=0,/silent)
	
		restore, ref_file[iref]

    if (qmask eq 'n') then mask = 1.

		ref_target = target
		ref_im1d = im1d
		ref_naxis3 = naxis3
	
		;is this star flagged to be bad, i.e. flag=0?
		ref_flag = flag[where(ref_target eq id)]

 		if (ref_target ne sci_target and ref_flag eq 1) then begin

;https://www.harrisgeospatial.com/docs/IDLmlafReLU.html
;https://pypi.org/project/pytorch-msssim/

			for i=0,sci_naxis3-1 do begin
;			for i=119,119 do begin
				for j=0,ref_naxis3-1 do begin

                    ssim[i,iref,j] = ms_ssim(sci_im[*,*,i]*mask, (ref_im[*,*,j]*mask))
                    if (finite(ssim[i,iref,j]) ne 1) then ssim[i,iref,j] = 0
                    ssim_file[i,iref,j] = ref_fitsfile[iref]
                    ssim_file_frame[i,iref,j] = j

				endfor

			endfor

		endif

		proceeding_text, loop=nfiles, i=iref, prompt='> Reference        '+string(iref+1,form='(I4)')

	endfor

	;---------------------------------------------------------------------
	
; 	if keyword_set(pca) then begin
	
    ;find the n-best reference frames for each science frame
    
    for i=0,sci_naxis3-1 do begin  ;loop over science frames

        for j=0,nbest-1 do begin ;loop to find XX best frames
        
            tmp = reform(ssim[i,*,*])
            dum = max(tmp,idx)
            ;print, dum, idx
            tmpidx = array_indices(tmp, idx)
            bestidx[i,j,*] = tmpidx
            best_ssim[i,j] = ssim[i,tmpidx[0],tmpidx[1]]
            ssim[i,tmpidx[0],tmpidx[1]] = -99
            best_ssim_file[i,j] = ssim_file[i,tmpidx[0],tmpidx[1]]
            best_ssim_file_frame[i,j] = ssim_file_frame[i,tmpidx[0],tmpidx[1]]

        endfor

    endfor
    
    st = {best_ssim_file:best_ssim_file, best_ssim_file_frame:best_ssim_file_frame, best_ssim:best_ssim}
    mwrfits, st, resdir+'Reference_'+file_id+'.fits', /create

endfor

print, ''
print, 'Done.'
print, ''
stop
end
