Usage of the RDI pipeline

Prerequisits
============

-have all final products (cosmetically reduced and image registered cubes and parallactic angles as txt file) of the NACO reduction in one base directory, e.g. NACO/<star1>/cube.fits, NACO/<star2>/cube.fits

-ISPY_RDI_prepare_data.pro
  -line 6: change path to the location of the individual observations (see previous item)
  -line 9: adjust results directory
  -line 12: define radial size of the images to be processed by RDI. The larger the slower the computation.
  -line 29 and 30: adjust the naming of the files, e.g. all image cubes have "....master.fits" as a common file name
  -line 167: definition of the radial size of the central mask in pixel, can be set to 0
  -line 168: definition of the mask applied to the edges of the image in pixel, can be set to 0
  -the resulting fits and sav files should be moved to a directory named "Lp/" if not already set in line 9

-prepare/modify text file "ISPY_RDI_targets.dat"
  -three columns per target: 1st col: exact target ID as defined by the object keyword in the fits header, 2nd col: type of target, e.g. PPD but is not mandatory and can be any non-empty string, 3rd col: flag that defines if the target can act as potential reference: 0=yes (can be reference), 1=no (cannot be reference). Targets with knwon CCs or beight disks are set to 1.
  
Finding the reference
=====================

-ISPY_RDI_find_ref
  -line 83: modify path. It is the directory where the output of ISPY_RDI_prepare_data.pro is stored
  -line 84: provide path to the location of the ISPY_RDI_targets.dat file
  -if you go through dozens of data sets it speeds things up if several parallel IDL sessions are started. Then you cal this routine and modify the range of the for loop in line 127. E.g. 1st instance: "for isci=0,9 do begin", 2nd instance: "for isci=10,19 do begin", ....
  -output is stored in a new directory "Reference_Lp/". It contains for each observation a fits file, which contains in the 1st extension the best SSIM value with the location of the corresponding reference file and the corresponding frame. Everything is sorted in descending order w.r.t. the SSIM value
  
-Subtraction of the reference PSF from the science frame using ISPY_RDI_PCAsubtraction.pro
  -line 5: modify the path to the directory containing the files after ISPY_RDI_prepare_data.pro was executed (should be "Lp/")
  -results are written into a new directory called "Results_Lp_PCA_nbest100/" and will contain for each observation follwoing files:
    -<id>_rdi_pca_median.fits final median combined image cube. 1st frame means 1 reference frame was used to construct a master reference, 2nd frame means 2 ref. frames were used, ....
    -<id>_rdi_pca_mean.fits mean combined final image
    -<id>_rdi_pca.fits  derotated ref. image subtracted cube before median combination
    -<id>_rdi_pca_median_convolved.fits same as <id>_rdi_pca_median.fits but images convolved with gaussian
    -<id>_rdi_pca_pca_median_<xx>RandomFrames_<yy>Modes.fits <xx> random frames of <yy> modes (= number of reference images used) to create median combined image, it is a simple check if signal remains if only a subset of images were used 
    -<id>_rdi_pca_median_<x>Frames_<yy>Modes.fits here only every 2nd, 3rd, 4th, 5th frame were used to create a median combined cube for a specific number of modes, it is a simple check if signal remains if only a subset of images were used 
    -<id>_rdi_pca_median_<1st, 2nd>HalfFrames_<yy>Modes.fits median combined image for the 1st or 2nd half of the observations, it is a simple check if signal remains if only a subset of images were used 
  
# RDI
