#!/usr/bin/python
# Filename: LOCI.py
import numpy as np

def MakePolarCoordinates(TargetImage): 
    npix = len(TargetImage)
    x = np.linspace(-npix/2, npix/2-1, npix)
    y = np.linspace(-npix/2, npix/2-1, npix)
    xv, yv = np.meshgrid(x, y)
  
    Rhos  = np.sqrt(xv**2 + yv**2)
    Thetas= np.arctan2(yv,xv)
    Thetas[Thetas<0]+=2*np.pi
    Thetas= np.nan_to_num(Thetas)
    
    return Rhos, np.transpose(Thetas)



# THIS IS A LITTERAL TRANSLATION FROM REMI'S MATHEMATICA CODE
# WITH SOME OPTIMIZATION (RUNS MUCH FASTER):
def LOCI(Target, Reference, rin, dr, Deltar, Thetamin, DeltaPhi, radialShiftBetweenSandO, Rhos, Thetas, RegNum, CriterionNumber, Mu):
    sizeImg    = len(Target)
    RhosFlat   = np.ndarray.flatten(Rhos)
    ThetasFlat = np.ndarray.flatten(Thetas)
    OzoneShape = np.zeros( (len(Target),len(Target) ))
    DropOnesList = []
    ChiSquareValues = []
    ReferenceFlat = []

    # FIRST, FIND OVERLAPPING REGIONS IN RHO AND THETA:
    
    #THIS IS MUCH MUCH FASTER THAN BEFORE (I.E. ELEMENT-BY-ELEMENT):
    idx1 = np.where( (Rhos>=rin) & (Rhos<(rin+dr)) ) 
    idx2 = np.where( (Thetas>=Thetamin) & (Thetas<(Thetamin + DeltaPhi)) )
    img=np.zeros( (sizeImg,sizeImg) )
    img[idx1]=1
    img[idx2]+=1
    Scoords = np.argwhere(img==2)

    #THIS IS MUCH MUCH FASTER THAN BEFORE (I.E. ELEMENT-BY-ELEMENT):
    idx1 = np.where( (Rhos>=rin + radialShiftBetweenSandO) & (Rhos<rin + radialShiftBetweenSandO + Deltar))
    idx2 = np.where( (Thetas>=Thetamin) & (Thetas< Thetamin + DeltaPhi)) 
    img=np.zeros( (sizeImg,sizeImg) )
    img[idx1]=1
    img[idx2]+=1
    Ocoords = np.argwhere(img==2)
    
    
    for i in xrange(len(Ocoords)):
        OzoneShape[Ocoords[i,0],Ocoords[i,1]]=1           
    Ocoords = np.transpose(np.where(OzoneShape==1))
    OT = Target[Ocoords[:,0],Ocoords[:,1]]
    ST = Target[Scoords[:,0],Scoords[:,1]]  

    
    for k in xrange(len(Reference)):
        OR = Reference[k][Ocoords[:,0],Ocoords[:,1]]
        ChiSquared = np.sum((OT - OR)**2)
        ChiSquareValues.append(ChiSquared)
        
    for k in xrange(len(Reference)):        
        if ChiSquareValues[k] > CriterionNumber*np.median(ChiSquareValues):
            DropOnesList.append(k)
    ReferenceFinal = np.delete(Reference, DropOnesList, axis=0)

    
    for k in xrange(len(ReferenceFinal)):
        ReferenceFlat.append(  np.ndarray.flatten(ReferenceFinal[k])  )    
    ReferenceFlat = np.transpose( np.asarray(ReferenceFlat) )
    
    #THIS IS MUCH MUCH FASTER THAN BEFORE (I.E. ELEMENT-BY-ELEMENT):
    idx1= np.where( (RhosFlat>=rin) & (RhosFlat<rin+dr))   
    idx2= np.where( (ThetasFlat>=Thetamin) & (ThetasFlat< Thetamin + DeltaPhi))
    img=np.zeros( (sizeImg*sizeImg) )
    img[idx1]=1
    img[idx2]+=1
    ScoordsFlat = np.argwhere(img==2)
    OcoordsFlat = ( np.transpose(np.where(np.ndarray.flatten(OzoneShape==1))) )

    # SECOND, THIS IS THE ACTUAL LOCI MATRIX ALGORITHM:
    OMatrix = ReferenceFlat[OcoordsFlat[:,0]]
    SMatrix = ReferenceFlat[ScoordsFlat[:,0]]
    
    AMatrix = np.dot(np.transpose(OMatrix),OMatrix) - Mu*np.dot( np.transpose(SMatrix),SMatrix)
    bvector = np.ndarray.flatten( np.dot( np.transpose(OMatrix),OT) - Mu*np.dot( np.transpose(SMatrix),ST) )
    
    ConditionNumber = np.ndarray.max(np.real(np.linalg.eigvals(AMatrix) ) ) / np.ndarray.min( np.real( np.linalg.eigvals(AMatrix) ) )
    try: Cks = np.dot( np.linalg.inv(AMatrix), bvector )
    except: 
        Cks = np.dot( np.linalg.inv( AMatrix + np.ndarray.max(np.real(np.linalg.eigvals(AMatrix)))*np.identity(len(AMatrix)) ), bvector)

    SR = np.dot( SMatrix,Cks )
    SI = np.transpose(ST - SR)
    
    return SI, Scoords



# THIS IS A LITTERAL TRANSLATION FROM REMI'S MATHEMATICA CODE:
def ApplyLOCI(Target,Reference, g, W, NA, dr,rstart, rend, radialShiftBetweenSandO, RegNum,CriterionNumber, Mu): 
    LengthOfImage =  len(Target)
    Area = NA*np.pi*(W/2.)**2
    DeltaR = np.sqrt(Area*g)
   
    FinalImage = np.zeros( (LengthOfImage, LengthOfImage) )
       
    Rhos, Thetas = MakePolarCoordinates(Target)
    radialbin = np.arange(rstart,rend,dr)

    for rin in radialbin:
        NAC = np.round( 2*np.pi/( g/2. + 2.*rin/W*(g/(np.pi*NA))**0.5 )**-1 ) #Number of Azimuthal Cuts
        DeltaPhi = 2*np.pi/NAC
        thetabin=np.arange(0,2*np.pi-DeltaPhi+DeltaPhi/100,DeltaPhi)
                
        for Thetamin in thetabin:
            SI, Scoords = LOCI(Target, Reference, rin, dr, DeltaR, Thetamin, DeltaPhi, radialShiftBetweenSandO, Rhos, Thetas, RegNum, CriterionNumber, Mu)
            for k in xrange(len(Scoords)):
                FinalImage[Scoords[k,0],Scoords[k,1]]=SI[k]
                     
    return FinalImage

