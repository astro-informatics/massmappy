from astropy.io import fits
import numpy as np
import scipy.io as sio
import os.path

# Download data from here https://des.ncsa.illinois.edu/releases/sva1/content
# You will need sva1_gold_r1.0_wlinfo.fits.gz, sva1_gold_r1.1_im3shape.fits.gz and sva1_gold_r1.0_bpz_point.fits.gz
# We would like to thank Chihway Chang and Joe Zunts for their help with constructing tha catalogue

info_check     = os.path.exists('data/sva1_gold_r1.0_wlinfo.fits.gz')
im3shape_check = os.path.exists('data/sva1_gold_r1.1_im3shape.fits.gz')
zinfo_check    = os.path.exists('data/sva1_gold_r1.0_bpz_point.fits.gz')

if not info_check:
	print "Cannot find data/sva1_gold_r1.0_wlinfo.fits.gz"
elif not im3shape_check:
	print "Cannot find data/sva1_gold_r1.1_im3shape.fits.gz"
elif not zinfo_check:
	print "Cannot find data/sva1_gold_r1.0_bpz_point.fits.gz"
else:

	print "Making the catalogue"
	print "This may take some time..."



	#read the info file
	cat = fits.open('data/sva1_gold_r1.0_wlinfo.fits.gz')


	ID       =  cat[1].data['COADD_OBJECTS_ID']
	RA       =  cat[1].data['RA']
	dec      =  cat[1].data['DEC']
	svflag   =  cat[1].data['SVA1_FLAG']
	flag     =  cat[1].data['IM3SHAPE_FLAG'] #im3shape flag, 0==ok 

	# #read the ellipticity catalogue
	cat = fits.open('data/sva1_gold_r1.1_im3shape.fits.gz')

	e1       =  cat[1].data['E_1']
	e2       =  cat[1].data['E_2']
	mcorr    =  cat[1].data['NBC_M']
	c1       =  cat[1].data['NBC_C1']
	c2       =  cat[1].data['NBC_C2']
	weight   =  cat[1].data['W']
	imerror  =  cat[1].data['ERROR_FLAG']
	imflag   =  cat[1].data['INFO_FLAG']

	# #now get redshifts 
	cat = fits.open('data/sva1_gold_r1.0_bpz_point.fits.gz')

	redshift =  cat[1].data['Z_MEAN']

	#now make selections based on table I here http://arxiv.org/pdf/1505.01871v2.pdf
	#note that the "conservative additive". In this paper https://arxiv.org/pdf/1504.03002.pdf they say: 
	#"we adopt the 'conservative additive' selection;
	#this results in small additive systematic uncertainties, but
	#possibly some moderate multiplicative systematic uncertainties.
	#For ngmix, this selection removes galaxies with S/N<20
	#and very small galaxies (Gaussian sigma smaller than the
	#pixel scale). For im3shape, it removes galaxies with S/N<15.
	#In both cases, there were many other selections applied to both
	#catalogs to remove stars, spurious detections, poor measurements,
	#and various other effects that significantly biased shear
	#estimates for both catalogs."


	good = np.logical_and(flag==0, dec>-70)
	good = np.logical_and(good, svflag==0)
	good = np.logical_and(good, dec<-40)
	good = np.logical_and(good, RA>60)
	good = np.logical_and(good, RA<95)
	good = np.logical_and(good, imerror==0.)
	good = np.logical_and(good, imflag==0.)
	good = np.logical_and(good, redshift>0.6)
	good = np.logical_and(good, redshift<1.2)

	RA     =  RA[good]
	dec    =  dec[good]
	e1     =  e1[good]
	e2     =  e2[good]
	mcorr  =  mcorr[good]
	c1     =  c1[good]
	c2     =  c2[good]
	weight =  weight[good]

	#number of objects
	print "Number of galaxies : ",  RA.size, "(it should be 793743)"

	sio.savemat('data/DES_test.mat',{'RA' : RA, 'dec' : dec, 'e1' : e1, 'e2' : e2, 'weight' : weight , 'mcorr' : mcorr , 'c1' : c1 , 'c2' : c2})
