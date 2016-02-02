#!/usr/bin/python
# Filename: AlignImages.py

import numpy as np
import scipy.interpolate as scint
import scipy.optimize as scopt

def interpolate(c, dim, tar, interp1, interp2):
    x = np.arange(0,dim)
    y = np.arange(0,dim)
    ref=np.zeros( [dim,dim] )
    mas=np.zeros( [dim,dim] )
 
    # Much faster than previous double FOR loop:
    mas = interp2(y-c[0], x-c[1])
    ref = c[2]*interp1(y-c[0], x-c[1])

    return np.log10( np.sum( ((ref-tar)*mas)**2 ) )


def align_images(dim, radius, reference, target):
    mask = np.ones( [dim,dim] )
    x    = np.arange(0,dim)
    y    = np.arange(0,dim)
    xmid = dim/2 
    ymid = dim/2     
    
    #binary mask centered on (xmid, ymid): 0 within "radius", 1 outside
    mask = np.array([ [0. if (np.sqrt( (i-xmid)**2 + (j-ymid)**2) < radius) else 1. for i in xrange(dim)] for j in xrange(dim)])

    alpha = 0.0
    beta  = 0.0
    nu    = 1.00
    guess = [alpha,beta,nu]

    interp_ref  = scint.RectBivariateSpline(x, y, reference, kx=3, ky=3)
    interp_mask = scint.RectBivariateSpline(x, y, mask, kx=1, ky=1)
   
    results = scopt.minimize(interpolate, guess, args=(dim,target,interp_ref,interp_mask), tol=1e-5)
    print "   alpha=%8.6f  beta=%8.6f  nu= %8.6f" %(results.x[0],results.x[1],results.x[2])

    # Much faster than previous double FOR loop:
    ref_aligned = results.x[2] * interp_ref(y-results.x[0], x-results.x[1])

    return ref_aligned


