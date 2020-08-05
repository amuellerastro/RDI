pro ISPY_RDI_PCAsubtraction

filter = 'Lp'

datapath = '/data/beegfs/astro-storage/groups/launhardt/amueller/NACO/DACE/data/'+filter+'/'

qmask = ''
read, 'Mask (y/n): ', qmask

nbest = 100
modes = dindgen(nbest)+1

if (qmask eq 'y') then refpath = datapath+'../Reference_'+filter+'_mask/' else refpath = datapath+'../Reference_'+filter+'_nomask/'

if (qmask eq 'y') then resdir = datapath+'../Results_'+filter+'_PCA_nbest'+strcompress(nbest,/rem)+'_mask/' else resdir = datapath+'../Results_'+filter+'_PCA_nbest'+strcompress(nbest,/rem)+'_nomask/'
spawn, 'mkdir -p '+resdir

nmodes = n_elements(modes)

;=========================================================================

if (filter eq 'Lp') then lambda = 3.8d-6

diam = 8.2
plsc = 27.19d-3
fwhm = lambda/diam*206265.d0/plsc

;=========================================================================

; file = file_search(datapath+'*fits', count=n)
reffile = file_search(refpath+'*fits', count=nref)


for xx=0,nref-1 do begin

    ;get file name
    fn = strmid(reffile[xx], strpos(reffile[xx],'Reference_',/reverse_search)+10, strlen(reffile[xx])-strpos(reffile[xx],'Reference_',/reverse_search)-15)
    
    file = datapath+fn+'.fits'
    print, file

    ;read in science
    im = mrdfits(file, 0, /silent)
    nframes = n_elements(im[0,0,*])
    dim = n_elements(im[*,0,0])
    
    paral = mrdfits(file, 1, /silent)
    
    diff_im = dblarr(dim, dim, nframes, nmodes)
    rot_diff_im = dblarr(dim, dim, nframes, nmodes)
    
    ;load reference
;     fn = strmid(file[xx], strpos(file[xx],'/',/reverse_search)+1, strlen(file[xx])-strpos(file[xx],'/',/reverse_search)-1)
;     fn = strmid(fn, 0, strpos(fn, '.', /reverse_search))
    st = mrdfits(reffile[xx], 1, /silent)
    best_ssim_file = st.best_ssim_file
    best_ssim_file_frame = st.best_ssim_file_frame
    best_ssim = st.best_ssim
    ;PCA subtraction of each science frame
    for i=0,nframes-1 do begin

        refframes = dblarr(dim, dim, nbest)
        for j=0,nbest-1 do begin

            ;workaround in case no reference files present, happend e.g. for Reference_CD-65149_2018-09-03T09_14_59_Lp.fits index 189,*
            ;using the references from previous frame at the moment
            ;if (strcompress(best_ssim_file[i,j],/rem) eq '') then begin

              ;tmp = mrdfits(strcompress(best_ssim_file[i-1,j],/rem),0,/silent)
              ;refframes[*,*,j] = tmp[*,*,best_ssim_file_frame[i-1,j]]
            
              ;if (best_ssim_file_frame[i-1,j] eq '') then begin
              
                ;tmp = mrdfits(strcompress(best_ssim_file[i-2,j],/rem),0,/silent)
                ;refframes[*,*,j] = tmp[*,*,best_ssim_file_frame[i-2,j]]
              
              ;endif

            ;endif else begin
        
              ;load reference frames
              tmp = mrdfits(strcompress(best_ssim_file[i,j],/rem),0,/silent)
              refframes[*,*,j] = tmp[*,*,best_ssim_file_frame[i,j]]
            
            ;endelse
    
        endfor
    
        ;PCA
        obj = refframes
        data = transpose(reform(obj,dim*dim,nbest))
        covMatrix = matrix_multiply(data, data, /btranspose)
        for k=0,nmodes-1 do begin
        
            eigenval = la_eigenql(covMatrix, EIGENVECTORS=eigenvect, range=[nbest-modes[k],nbest-1], /double)
            eigenval = reverse(eigenval)
            if ((size(eigenvect))[0] gt 1) then eigenvect = reverse(eigenvect,2) else eigenvect = reverse(eigenvect)
            pc_orig = matrix_multiply(eigenvect,data,/atranspose)
            pc = pc_orig
            
            for l=0,modes[k]-1 do pc[l,*] = pc_orig[l,*]/(eigenval[l])
        
            scidata = transpose(reform(im[*,*,i],dim*dim))
            s1 = matrix_multiply(pc_orig, scidata, /btranspose)
            ref_frame = reform(matrix_multiply(s1, pc), dim, dim)
            diff_im[*,*,i,k] = im[*,*,i]-ref_frame
            rot_diff_im[*,*,i,k] = rot(reform(diff_im[*,*,i,k]), -paral[i], 1.0, dim/2., dim/2., cubic=-0.5, /pivot)
            
        endfor
        
    endfor
    
    medim = dblarr(dim,dim,nmodes)
    meanim = dblarr(dim,dim,nmodes)
    cimg = dblarr(dim,dim,nmodes)
    for i=0,nmodes-1 do begin
    
        medim[*,*,i] = median(rot_diff_im[*,*,*,i], dim=3)
        meanim[*,*,i] = mean(rot_diff_im[*,*,*,i], dim=3)
        cimg[*,*,i] = filter_image(medim[*,*,i], fwhm_gaussian=0.5*fwhm)
        
    endfor
    
    writefits, resdir+fn+'_rdi_pca_median.fits', medim
    writefits, resdir+fn+'_rdi_pca_mean.fits', meanim
    writefits, resdir+fn+'_rdi_pca.fits', rot_diff_im
    writefits, resdir+fn+'_rdi_pca_median_convolved.fits', cimg
    
    ;=====================================================================

    ;Subsets
    submode = round(nmodes/2)

    ;choose every n-th frame
    nth = [1,2,3,4,5]
    for i=0,n_elements(nth)-1 do begin

	    if (nframes gt nth[i]) then begin

            idx = round(linspace(0,nframes-1,nframes/double(nth[i])))
            
            if (n_elements(idx) gt 1 and idx[0] ne -1) then begin
                    
                newim = median(rot_diff_im[*,*,idx,submode-1], dim=3, /even)
                writefits, resdir+fn+'_rdi_pca_median_'+strcompress(nth[i],/rem)+'Frames_'+strcompress(submode,/rem)+'Modes.fits', newim, hdr

            endif

	    endif

    endfor
    
    ;first half and second half
    if (nframes ge 8.) then begin
    
        writefits, resdir+fn+'_rdi_pca_median_1stHalfFrames_'+strcompress(submode,/rem)+'Modes.fits', median(rot_diff_im[*,*,0:round(nframes/2.)-1, submode-1], dim=3, /even), hdr
        writefits, resdir+fn+'_rdi_pca_median_2ndHalfFrames_'+strcompress(submode,/rem)+'Modes.fits', median(rot_diff_im[*,*,round(nframes/2.):*,submode-1], dim=3, /even), hdr

        ;choose randomly 20%, 35%, ... of the frames and do this niter times
        niter = 100
        frac = [0.2,0.35,0.5,0.75]
        for i=0,n_elements(frac)-1 do begin
        
            imr = dblarr(n_elements(rot_diff_im[*,0,0,0]),n_elements(rot_diff_im[0,*,0,0]),niter)
    
            for j=0,niter-1 do begin
    
                n = round(frac[i]*nframes)
                tmp = round(nframes*randomu(seed,n))
                idxs = tmp[sort(tmp)]
                if (n_elements(idxs) gt 1 and idx[0] ne -1) then imr[*,*,j] = median(rot_diff_im[*,*,idxs,submode-1], dim=3, /even)

            endfor
            
            if (n_elements(idxs) gt 1 and idx[0] ne -1) then writefits, resdir+fn+'_rdi_pca_median_'+sigfig(frac[i],2)+'RandomFrames_'+strcompress(submode,/rem)+'Modes.fits', imr, hdr

        endfor
    
    endif
    
    proceeding_text, loop=nref, i=xx, prompt='> Target        '+string(xx+1,form='(I4)')
    
endfor

stop
end
