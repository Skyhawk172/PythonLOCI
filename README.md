# PythonLOCI


## alignImage.ipynb
This ipython notebook was used for development of the LOCI algorithm in
python. I got it to a point where everything matched very nicely the results
from Remi's Mathematica notebook, so I decided to create modules out of it so
that I can use it from anywhere. See below:



##ApplyLOCI.py
Script that runs the AlignImages and LOCI modules on the specified data. 
The script uses command-line arguments and, by design, expects the PSF files to
contain a string with "runXX" in it in order to process a given run number. By
default, it will process all the runs in the specified directory, unless the
--run flag is used.

To run it, simply type: python ApplyLOCI.py directory_to_PSFs [--run RunNumber]

For more information: python ApplyLOCI.py --help



##LOCI.py
Actual LOCI algorithm. This runs much faster than the Mathematica notebook.
This function is called from within ApplyLOCI.py and takes for arguments:
 * input target image
 * input reference images (typically, 4, 9, or 25)
 * kwargs for specifying various parameters of the LOCI algorithm.



##AlignImages.py
Cubic least-square minimization technique to align the reference images with the
target image. This version also runs much faster than before.
This function is called from within ApplyLOCI.py and takes for arguments:
 * sizeImg: the size of the image (integer)
 * radius: the radius of the region to NOT use for the minimization
 * target: the image of the science target
 * refs[i]: reference image to align with the target
