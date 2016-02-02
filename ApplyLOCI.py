import alignImages
import LOCI
import numpy as np
import pyfits
import argparse
import os, glob


def makeHeader(hdr, directory, target, references, LOCIargs ):
    # To keep track in FITS header of what files were used for LOCI    
    hdr.add_history(os.getcwd())
    hdr.add_history(target[0])
    for iref in references:
        hdr.add_history(iref)
        
    for key, value in LOCIargs.iteritems():
        hdr.add_history(key+"="+str(value))

    return hdr


def process_run(nruns, **kwargs):
    info= '\nStarting LOCI run %d -->' %nruns
    print info,

    targetfile=glob.glob("*run"+str(nruns)+"_"+"*ScienceTarget*"+"*.fits")

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
        refs_aligned[i]=alignImages.align_images(sizeImg,radius,refs[i],target)

    unoccultedfile=glob.glob("*run"+str(nruns)+"_"+"*Unocculted*.fits")
    hdu=pyfits.open(unoccultedfile[0])
    unocculted=hdu[0].data

    finalImage=LOCI.ApplyLOCI(target, refs_aligned, **kwargs)
    finalImage=finalImage/np.max(unocculted)

    hdu = pyfits.PrimaryHDU(finalImage)
    hdr = hdu.header
    header = makeHeader(hdr, directory, targetfile, reffile, LOCIargs)
    outname = "Map_LOCI_run%d.fits" %nruns
    hdu.writeto(outname,clobber=True)
    print "Output file:", outname

    return




#####################################################################
# MAIN Method
#####################################################################
if __name__ == "__main__":

    #-------------------------------
    # COMMAND-LINE ARGUMENT PARSER:
    #-------------------------------
    parser = argparse.ArgumentParser(description="Apply LOCI to sets of PSFs")

    parser.add_argument("dir"   , type=str, help="Input directory (absolute or relative)")
    parser.add_argument("--run" , default="all", help="run number (default = all)")

    args = parser.parse_args()
    directory = args.dir


    #----------------------------
    # DEFAULT IMAGE PARAMETERS:
    #----------------------------
    if "MIRI" in directory:
        sizeImg  = 64
        radius   = 25
        pixelSize= 0.11 #arcseconds
    elif "NIRCam" in directory:
        if "F210" in directory:
            print "NIRCam Short Wavelength"
            sizeImg  = 222 
            radius   = 70 
            pixelSize= 0.032 #0.065 #arcseconds
        else:
            print "NIRCam Long Wavelength"
            sizeImg  = 109 
            radius   = 25 
            pixelSize= 0.065 #arcseconds


    #----------------------------
    #ARGUMETS FOR LOCI ALGORITHM:
    #----------------------------
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


    #----------------------------
    # RUN LOCI ON EACH RUN:
    #----------------------------
    os.chdir(directory)
    if args.run == "all":
        nruns = len(glob.glob("*ScienceTarget*.fits"))
        for irun in xrange(1, nruns+1):
            process_run(irun, **LOCIargs)

    else:
        nruns = int(args.run)
        process_run(nruns, **LOCIargs)
        



