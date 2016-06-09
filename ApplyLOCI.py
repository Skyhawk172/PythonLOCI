#################################################################
# WHAT: -script to apply LOCI on a set of reference images and a
#        science target. Also applies a classical subtraction
#        between the science target and another reference target.
#
#       - requires the LOCI.py and AlignImages.py modules (see
#         Github/Skyhawk172).
#
# HOW: python ApplyLOCI.py --help
#
# WHO: C-P Lajoie - STScI
#
# WHEN: February 2016
#################################################################

import AlignImages
import LOCI
import numpy as np
import pyfits
import argparse
import os, glob, sys
import matplotlib.pyplot as plt


def makeHeader(hdr, directory, target, references, **kwargs):
    # To keep track in FITS header of what files were used for LOCI    
    hdr.add_history(os.getcwd())
    hdr.add_history(target[0])
    for iref in references:
        hdr.add_history(iref)
        
    for key, value in kwargs.iteritems():
        hdr.add_history(key+"="+str(value))

    return hdr


def process_run(nruns, mask, **kwargs):
    info= '\nStarting LOCI run %02d -->' %nruns
    print info

    targetfile= glob.glob("*run"+str(nruns)+"_"+"*ScienceTarget*"+"*.fits")


    #------------------------------
    # APPLY LOCI FIRST:
    #------------------------------
    hdu=pyfits.open(targetfile[0])
    target=hdu[0].data

    reffile=glob.glob("*run"+str(nruns)+"_"+"*dither*.fits")
    nimages=len(reffile)

    refs=np.zeros( (nimages,sizeImg,sizeImg) )
    refs_aligned=np.zeros( (nimages,sizeImg,sizeImg) )
    for i in xrange(nimages):
        print i+1, reffile[i]
        hdu=pyfits.open(reffile[i])
        refs[i]=hdu[0].data
        refs_aligned[i]=AlignImages.align_images(sizeImg, mask, refs[i], target)

    #unoccultedfile=glob.glob("*run"+str(nruns)+"_"+"*Unocculted*.fits")
    #hdu=pyfits.open(unoccultedfile[0])
    #unocculted=hdu[0].data

    finalImage=LOCI.ApplyLOCI(target, refs_aligned, **kwargs)

    hdu = pyfits.PrimaryHDU(finalImage)
    hdr = hdu.header
    header = makeHeader(hdr, directory, targetfile, reffile, **kwargs)#LOCIargs)
    outname = "Map_LOCI_run%02d.fits" %nruns
    hdu.writeto(outname,clobber=True)
    print "Output file:", outname



    #------------------------------
    # APPLY CLASSICAL SUBTRACTION:
    #------------------------------
    reffile = glob.glob("*run"+str(nruns)+"_"+"*ReferenceTarget*"+"*.fits")
    hdu=pyfits.open(reffile[0])
    reference=hdu[0].data

    refaligned=AlignImages.align_images(sizeImg, mask, reference, target)

    finalImage=target-refaligned

    hdu = pyfits.PrimaryHDU(finalImage)
    hdr = hdu.header
    header = makeHeader(hdr, directory, targetfile, reffile, **kwargs)#LOCIargs)
    outname = "Map_CLAS_run%02d.fits" %nruns
    hdu.writeto(outname,clobber=True)
    print "Output file:", outname



    return




#-----------------------------------#
# MAIN Method
#-----------------------------------#
if __name__ == "__main__":

    #-------------------------------#
    # COMMAND-LINE ARGUMENT PARSER:
    #-------------------------------#
    parser = argparse.ArgumentParser(description="Apply LOCI to sets of PSFs")

    parser.add_argument("dir"      , type=str, help="Input directory (absolute or relative)")
    parser.add_argument("--run"    , default="all", help="run number (default = all)")
    parser.add_argument("--imgSize", type=int, help="image size (pixels)")
    parser.add_argument("--rad"    , type=int, help="radius to ignore for alignment function (pixels)")
    parser.add_argument("--pixSize", type=int, help="instrument's pixel size (arcseconds)")
    args = parser.parse_args()
    directory = args.dir


    #-------------------------------#
    # DEFAULT IMAGE PARAMETERS:
    #-------------------------------#
    if "Results/MIRI" in directory:
        sizeImg  = 64   if args.imgSize==None else args.imgSize
        radius   = 25   if args.rad    ==None else args.rad
        pixelSize= 0.11 if args.pixSize==None else args.pixSize 
    elif "Results/NIRCam" in directory:
        if "F210M" in directory:
            print "\nNIRCam Short Wavelength"
            sizeImg  = 222   if args.imgSize==None else args.imgSize
            radius   = 75    if args.rad    ==None else args.rad
            pixelSize= 0.032 if args.pixSize==None else args.pixSize
        else:
            print "\nNIRCam Long Wavelength"
            sizeImg  = 109   if args.imgSize==None else args.imgSize
            radius   = 40    if args.rad    ==None else args.rad
            pixelSize= 0.065 if args.pixSize==None else args.pixSize


    # MAKE ALIGNMENT MASK FOR EACH CORON. CONFIGURATION:
    xmid = sizeImg/2 
    ymid = sizeImg/2   

    mask = np.array([ [0. if (np.sqrt( (i-xmid)**2 + (j-ymid)**2) < radius) else 1. for i in xrange(sizeImg)] for j in xrange(sizeImg)])
    if "MIRI" in directory: 
        mask[ymid-4:ymid+5, :] = 0.
        mask[:, xmid-4:xmid+5] = 0.
    if "MASKSWB" in directory: mask[ymid-8:ymid+9, :] = 0.
    if "MASKLWB" in directory: mask[ymid-8:ymid+9, :] = 0.

    plt.imshow(mask)
    raw_input()
    #-------------------------------#
    #ARGUMENTS FOR LOCI ALGORITHM:
    #-------------------------------#
    g  = 1.75
    NA = 400
    dr = 1.25
    CriterionNumber = np.inf
    Mu= 0.0 
    W = (10.*10**-6)/6.5 * 180./np.pi * 3600./pixelSize
    rstart = 0.5
    rend = sizeImg/2. #30 for MIRI originally
    RegNum = 10**-9
    radialShiftBetweenSandO = 0.0

    LOCIargs = {"g":g, "NA":NA, "dr":dr, "CriterionNumber":CriterionNumber, "Mu":Mu, "W":W, 
            "rstart":rstart, "rend":rend, "RegNum":RegNum, "radialShiftBetweenSandO":radialShiftBetweenSandO}



    #-------------------------------#
    # RUN LOCI ON EACH RUN:
    #-------------------------------#
    os.chdir(directory)
    if args.run == "all":
        nruns = len(glob.glob("*ScienceTarget*.fits"))
        for irun in xrange(1, nruns+1):
            process_run(irun, mask, **LOCIargs)

    else:
        nruns = int(args.run)
        process_run(nruns, mask, **LOCIargs)
        



