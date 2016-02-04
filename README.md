# PythonLOCI

My conversion of a Mathematica notebook to Python for the LOCI and alignment
algorithms.

The basic idea here is that the user provides (1) one science target image, (2)
numerous reference images (dithered; typically 9), and (3) an unocculted image
(i.e. no coronagraph mask, for normalization purposes).  The code aligns all the
reference images with the target image and then applies the LOCI algorithm to
the set of images to create a "synthetic" reference image that optimizes the PSF
subtraction.	  


## alignImage.ipynb 

This ipython notebook was used for development of the LOCI algorithm in
python. I got it to a point where everything matched very nicely the results
from Remi's Mathematica notebook, so I decided to create modules out of it so
that I can use it from anywhere. See below.



##ApplyLOCI.py

This script runs the "AlignImages.py" and "LOCI.py" modules on the specified
data. By design, this is the only script you should have to interact with as the
other two are called from within it.

The script uses command-line arguments and, by design, expects the PSF
files to contain a string with "runXX" in it in order to process a given run
number. Moreover, the script expect the target image filename to contain the
string "ScienceTarget", the reference images to contain the string "dither", and
the unocculted image to contain the string "Unocculted". Here are a few example
filenames:
* PSF_1065_run7_ScienceTarget_204.fits 
* PSF_1140_run3_Unocculted_204.fits
* PSF_210_run1_Reference_dither1_136.fits

The user may need to customize these strings by simply searching for "glob"
statements.

By default, it will process all the runs in the specified directory,
unless the --run flag is used. The command-line arguments are (-- for optional):

* dir: full path or relative path to directory containing FITS files.
* --run: run number (default: all)
* --imgSize: size of images in pixels (default: 64) 
* --pixSize: instrument's pixel size in arcseconds (instrument dependent)
* --rad: radius from center of image to ignore during the alignment routine (instrument dependent).

Note also that this script is somewhat custom built, especially when it comes to
"finding" the files in the directory and/or determining the instrument/filter
being used. The user might have to modify the code accordingly. 

For more information: python ApplyLOCI.py --help



##LOCI.py

Actual LOCI algorithm. This runs much faster than the Mathematica notebook.
This function is called from within ApplyLOCI.py and takes for arguments:
 * input target image
 * list of input reference images (typically, 4, 9, or 25)
 * kwargs for specifying various parameters of the LOCI algorithm, which are 
   currently hardcoded (see ApplyLOCI.py)



##AlignImages.py

Cubic least-square minimization technique to align the reference images with the
target image. This version also runs much faster than before.
This function is called from within ApplyLOCI.py and takes for arguments:
 * sizeImg: the size of the image (integer)
 * radius: the radius of the region to NOT use for the minimization
 * target: the image of the science target
 * refs[i]: reference image to align with the target
