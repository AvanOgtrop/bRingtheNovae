################################################################################
#Code developed by Aniek van Ogtrop                                            #
#Last edit: 2020 September                                                     #
#Functions used in FullPipeline.py                                             #
################################################################################

#necessary imports
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import h5py
import datetime
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip

import glob
import pickle

import sys
sys.path.append('/disks/web1/bring/products/')
sys.path.append('/net/belterwijde/data2/stuik/code')
sys.path.append("/net/belterwijde/data2/vogtrop")

import burp
from bringreduce import astrometry
from bringreduce import photometry
from bringreduce import configuration as cfg

"""
This file contains the functions used in FullPipeline.py to extract a light 
curve of a star when only the RA and DEC are known.

Functions included in this file:
1. get_to_degrees(): change the RA and DEC from hour notation to degrees

2. find_3_closest(): find three closest reference stars that meet the criteria

3. name_get(): get the filename of the files with all three reference stars

4. get_lstseq(): find the lstseq where all three reference stars are present

5. conv(): convert the RA and DEC of the target to x and y coordinates

6. RA_DEC_conversion(): collect the x and y coordinates of the target star

7. cube(): make a post stamp selection of the images around the target star

8. find_flux(): find the flux of the target and reference stars

9. make_mask(): make a mask that excludes the period of nova eruption

10. periodicity(): filter out the periodicities provided

11. calibration(): apply the full calibration with all periodicities
"""

#to convert the RA and DEC to degrees
def get_to_degrees(RA_hr,RA_min, RA_sec,DEC_deg,DEC_arcmin,DEC_arcsec):
    """get_to_degrees() changes the RA and DEC from hour minute notation to 
degrees.
INPUT:
RA_hr: the right ascension hours
RA_min: the right ascension minutes
RA_sec: the right ascension seconds
DEC_deg: the declination degrees
DEC_arcmin: the declination arcminutes
DEC_arcsec: the declination arcseconds
OUTPUT:
RA: the right ascension in degrees
DEC: the declination in degrees
"""
    RA=(RA_hr+(RA_min/60.)+(RA_sec/3600.))*15
    DEC=abs(DEC_deg)+DEC_arcmin/60.+DEC_arcsec/3600.
    
    #to make sure the declination is negative if need be
    if DEC_deg<0:
        DEC*=-1
    return(RA,DEC)


#Find the three closest stars that meet the criteria
def find_3_closest(RA,DEC,starcatsel):
    """find_3_closest() finds the three closest stars to be used as reference 
stars. Stars that are variable or have a blend larger than 0.1 are rejected as
well as stars too close to the target or the other reference stars.
INPUT:
RA: the right ascension of the target star in degrees
DEC: the declination of the target star in degrees
starcatsel: the star catalogue to be used
OUTPUT:
idx: array containing the indices of the selected reference stars in starcatsel
ASCC: array containing the ASCC names of the reference stars (string)
HD: array containing the HD names of the reference stars (string)     
ra: array containing the RA of the target star and the reference stars
dec: array containing the DEC of the target star and the reference stars
Vmag: array containing the V magnitudes of the reference stars
    """
    #To find the closest 3 stars:
    c_target=SkyCoord(RA,DEC,unit="deg")
    c_catalogue=SkyCoord(starcatsel["_RAJ2000"],starcatsel["_DEJ2000"],unit="deg")
    deltaDist=c_target.separation(c_catalogue).deg
    
    #creating a new sorted array with the ASCC name in the first column and the 
    #deltaDist in the second column ascending
    sortedidx=np.argsort(deltaDist)
    sortedcat=starcatsel[sortedidx]
    
    #setting initial values and empty arrays to be filled later
    approved=0
    i=0
    
    idx=np.array([],dtype="int")
    ASCC=np.array([],dtype="int")
    HD=np.array([],dtype="int")
    ra=np.array([RA])
    dec=np.array([DEC])
    Vmag=np.array([])
    c_ref=np.array([])
    
    #loop over the stars (in order of closeness) until three stars are found 
    #that meet the requirements:
        #1. the star is not the target star/the star is not too close to the target
        #2. the star is not a variable star
        #3. the star has a blend of less than 0.1
        #4: the star is not too close to the other reference stars

    #if after 150 stars the three reference stars aren't found the loop is ended
    while approved<3 and i<150:
        
        if float(deltaDist[sortedidx][i])<1:
            print ("The star is too close to the target star")
        elif sortedcat[i]["Var"]!=0:
            print ("This star: {}, is variable".format(sortedcat[i]["ASCC"]))
        elif sortedcat[i]["Blend"]>0.1:
            print ("This star: {}, is blended".format(sortedcat[i]["ASCC"]))
        elif (len(c_ref)>0) and (np.min(c_ref.separation(SkyCoord(starcatsel["_RAJ2000"][sortedidx][i],starcatsel["_DEJ2000"][sortedidx][i],unit="deg")).deg)<0.5):
            print ("This star is too close to another reference star")
        else:
            print (sortedcat[i]["_RAJ2000"],sortedcat[i]["_DEJ2000"],sortedcat[i]["Vmag"])
            approved+=1
            idx=np.append(idx,sortedidx[int(i)]) 
            ASCC=np.append(ASCC,sortedcat[i]["ASCC"])
            HD=np.append(HD,sortedcat[i]["HD"])
            ra=np.append(ra,sortedcat[i]["_RAJ2000"])
            dec=np.append(dec,sortedcat[i]["_DEJ2000"])
            Vmag=np.append(Vmag,sortedcat[i]["Vmag"])
            c_ref=SkyCoord(ra,dec,unit="deg")
        i+=1
        
    return (idx,ASCC,HD,ra,dec,Vmag)


#get the file names containing all three reference stars
def name_get(curcam,curlst,x,y,binnedLUT):
    """name_get() checks if the camera exists and if the image exists and then
gives the name of the file the image is in and the lstseq of that observation. 
INPUT:
curcam: index of the camera as per burp.cams
curlst: the lstseq of the image (so even 50s)
x: list of the x coordinates of the reference stars
y: list of the y coordinates of the reference stars   
binnedLUT: binnedLUT pickle file
OUTPUT:
curlst+24: the lstseq of the observation for the long exposure
binnamelong: the name of the file corresponding to the lstseq of the long 
             exposure
    """

    cameras = ['SAW', 'SAE', 'AUW', 'AUE']
    
    #only continue when camera is there
    if not np.isin(burp.cams[curcam],cameras):
        print ("This camera is not available")
        return (0,0)
    
    #set the right camera
    cam=burp.cams[curcam]
    curcamidx=cameras.index(burp.cams[curcam])
    
    #from binnedLUT the date is extracted
    datestr = binnedLUT[curcamidx][np.uint64(curlst)]
    
    #if there is no file at that curlst, binnedLUT returns -1
    if int(datestr) == -1: 
        print('Binned image not found')
        
        return (datestr,0)
        
    
    #NOTE: this next while loop is only necessary if images are taken into 
    #account that don't have 50 images in the binned period
    else:
        midlstseq=24

        basepath = '/net/belterwijde/data2/bring/products/'
        listpaths=glob.glob(basepath+datestr+cam+'/binned/bin_'+"*.fits")
        
        filefound=False
        while filefound==False:
            
            binnamelong=basepath+datestr+cam+'/binned/bin_'+str(curlst+midlstseq)+cam+ '.fits'
            if np.isin(binnamelong,listpaths):
                filefound=True
                
                if (curlst+midlstseq)%2!=0:
                    binnamelong = basepath+datestr+cam+'/binned/bin_'+str(curlst+midlstseq-1)+cam+ '.fits'
                
            #this will only happen if the set is not made up of 50 images
            else:
                midlstseq-=1
      
        return (curlst+24,binnamelong)
   
        
#find all lstseq where the three reference stars are present        
def get_lstseq(dir,detrendStars,camera, star):
    """get_lstseq() compares the detrend data lstseq of the reference stars and
filters to only include the lstseq that are present in all reference star data. 
From this it finds the name of the file the image is in and saves this to a file
with the name "files_{star}_{camera}.txt".
INPUT:
dir: string containing the path to the working directory
detrendStars: an array with the detrend data of the reference stars
camera: index of the camera as per burp.cams
star: name of the target star
OUTPUT:
lstseq+24: an array containing the lstseq of the data where all reference stars
        are present
x: an array containing the x coordinates of the reference stars
y: an array containing the y coordinates of the reference stars
jd: an array containing the Julian date of the data where all reference stars 
    are present
lst: an array containing the lst of the data where all reference stars are 
     present
names: an array containing the file names where all reference stars are present    
    """
    
    #the catalogue that tells if the date has been observed
    catpath =  '/net/belterwijde/data2/bring/Catalogues/'
    binnedLUT = pickle.load(open(catpath+'/binnedLUT.pkl','r'))[0]
    
    #empty lists to be filled later
    names=[]
    Lstseq=[]
    X=[]
    Y=[]
    Jd=[]
    Lst=[]

    #only select the detrend data with the chosen camera and where 50 
    #observations are made
    for j, detrend in enumerate(detrendStars):
        mask=(detrend["camera"]==camera)*(detrend["nobs"]==50)
        Lstseq.append(detrend["lstseq"][mask])
        X.append(detrend["x"][mask])
        Y.append(detrend["y"][mask])
        Jd.append(detrend["jd"][mask])
        Lst.append(detrend["lst"][mask])

    #flatten list to find the unique values and the counts of each of them
    flat_lstseq= [item for sublist in Lstseq for item in sublist]
    unique_values, counts = np.unique(flat_lstseq, return_counts=True)
    
    #select the values that occur for each star
    mask=np.where(counts==len(detrendStars))[0]
    sel_lstseq=unique_values[mask]
    
    numbStars=len(detrendStars)
    numbLstseq=len(sel_lstseq)
    
    #create empty arrays in the right shape
    lstseq=np.array([None]*numbStars*numbLstseq).reshape(numbLstseq,numbStars)
    x=np.array([None]*numbStars*numbLstseq).reshape(numbLstseq,numbStars)
    y=np.array([None]*numbStars*numbLstseq).reshape(numbLstseq,numbStars)
    jd=np.array([None]*numbStars*numbLstseq).reshape(numbLstseq,numbStars)
    lst=np.array([None]*numbStars*numbLstseq).reshape(numbLstseq,numbStars)
    
    #fill the arrays with the selected values
    for i in range(len(detrendStars)):
        mask=np.isin(Lstseq[i],sel_lstseq)
        
        lstseq[:,i]=np.array(Lstseq[i])[mask]
        x[:,i]=np.array(X[i])[mask]
        y[:,i]=np.array(Y[i])[mask]
        jd[:,i]=np.array(Jd[i])[mask]
        lst[:,i]=np.array(Lst[i])[mask]
    
    #initiate lists to collect the lstseq
    t=[]
    
    #to clear the file
    open("{}/files_{}_{}.txt".format(dir,star,camera), "w").close()
    
    #collect the name of the file for each lstseq
    for i in range(numbLstseq): 
        xsel=x[i]
        ysel=y[i]
 
        T,binname=name_get(camera,sel_lstseq[i],xsel,ysel,binnedLUT)
        
        #only if it worked append the value
        if T>0:
            t.append(T)
        names.append(binname)
      
    names=np.array(names)
    mask=np.where(names!="0")[0]
    lstseq=lstseq[:,0][mask]
    x=x[mask]
    y=y[mask]
    jd=jd[mask]
    lst=lst[mask]
    names=names[mask]
    
    #save the information to files
    np.savetxt("{}/files_{}_{}.txt".format(dir,star,camera),names,fmt="%s")
    np.savetxt("{}/x_reference_{}_{}.txt".format(dir,star,camera),x,fmt="%f")
    np.savetxt("{}/y_reference_{}_{}.txt".format(dir,star,camera),y,fmt="%f")
    np.savetxt("{}/jd_{}_{}.txt".format(dir,star,camera),jd[:,0],fmt="%f")
    np.savetxt("{}/lstseq_{}_{}.txt".format(dir,star,camera),lstseq+24,fmt="%f") #changed
    return (lstseq+24,x,y,jd[:,0],lst[:,0],names)


#convert the RA and DEC to x and y for all files
def conv(ra, dec, curlstseq, curlst, jd, astrofile, fast, camera):
    """conversion() converts the given ra and dec to x and y positions on the 
files given.
INPUT:
ra: an array containing the right ascension of the target star(s)
dec: an array containing the right ascension of the target star(s)
curlstseq: the lstseq of the image
curlst: the lst of the image
jd: the Julian date of the image
astrofile: string of the path of the fast photometry file
fast: the fast photometry file
camera: index of the camera as per burp.cams
OUTPUT:
x0: the x coordinate on the file corresponding to the ra and dec given
y0: the y coordinate on the file corresponding to the ra and dec given
err0: gives False if the star is outside the image or too close to the edges,
      gives True if the star is in the image and the ra and dec are succesfully
      converted
    """
    
    nx = 4008
    ny = 2672
    
    curastro = np.where((fast['astrometry/lstseq'][()] / 50) == (curlstseq / 50))[0][0]
    order = fast['astrometry/x_wcs2pix'][curastro].shape[0]-1
    astrolst = fast['station/lst'][np.where(fast['station/lstseq'][()] == (fast['astrometry/lstseq'][curastro]))[0][0]]
    
    wcspars = {'crval' : fast['astrometry/crval'][curastro].copy(),
               'crpix' : fast['astrometry/crpix'][curastro].copy(),
               'cdelt' : [0.02148591731740587, 0.02148591731740587],
               'pc'    : fast['astrometry/pc'][curastro].copy(),
               'lst'   : astrolst}
               
    polpars = {'x_wcs2pix' : fast['astrometry/x_wcs2pix'][curastro].copy(),
               'y_wcs2pix' : fast['astrometry/y_wcs2pix'][curastro].copy(),
               'x_pix2wcs' : fast['astrometry/x_pix2wcs'][curastro].copy(),
               'y_pix2wcs' : fast['astrometry/y_pix2wcs'][curastro].copy(),
               'nx'    : nx,
               'ny'    : ny,
               'order' : order}
               
    astro = astrometry.Astrometry(wcspars, polpars)
    #due to an implementation error of the equinox for the cameras in
    #South Africa before 2018 a distinction is made to get the
    #correct conversion
    if (int(astrofile[-16:-12]) < 2018) & (camera[:2] == 'SA'):
        x0, y0, err0 = astro.world2pix(curlst, ra, dec)
    else:
        x0, y0, err0 = astro.world2pix(curlst, ra, dec, jd=jd)
    return (x0, y0, err0)

#collect the x and y of the target star       
def RA_DEC_conversion(dir,RA_target,DEC_target,files,star,camera,lstseq,lst,jd,xref,yref):
    """RA_DEC_conversion() gets the x and y coordinates for all images in files
from the RA and DEC. The x and y is saved in a file with the name 
"x_y_{star}_{camera}.txt" with the first column the lstseq and the next x and y.
INPUT:
dir: string containing the path to the working directory
RA_target: the right ascension of the target in degrees
DEC_target: the declination of the target in degrees
files: an array containing the paths of the images
star: name of the target star
camera: index of the camera as per burp.cams
lstseq: an array containing the lstseq
lst: an array containing the lst
jd: an array containing the Julian date
xref: an array containing the x coordinates of the reference stars
yref: an array containing the y coordinates of the reference stars
OUTPUT:
x_target: an array containing the x coordinates on the images
y_target: an array containing the y coordinates on the images
files: an array containing the paths of the images
lstseq: an array containing the lstseq
lst: an array containing the lst
jd: an array containing the jd
xref: an array containing the x coordinates of the reference stars
yref: an array containing the y coordinates of the reference stars               
    """
    
    #empty arrays to be filled later
    x_target=np.array([])
    y_target=np.array([])
    
    #emptying the file
    open("{}/x_y_{}_{}.txt".format(dir,star,camera), "w").close()
    
    camstr = ['AUE', 'AUW', 'LSC', 'LSE', 'LSN', 'LSS', 'LSW', 'SAE', 'SAW', 'LPC', 'LPE', 'LPN', 'LPS', 'LPW']
    
    #collecting all fast photometry files
    fastfiles = np.sort(glob.glob("/net/belterwijde/data3/bring/products/"+'*'+camstr[camera]+'/lightcurves/fast*'))
    
    
    #build an index
    minlstseq = np.zeros(len(fastfiles))
    maxlstseq = np.zeros(len(fastfiles))
    for fff, jjj in zip(fastfiles, range(len(fastfiles))):
        try:
            f = h5py.File(fff, 'r')
            try:
                station = f['station']
                if len(np.where(np.array(station.keys()) == 'LSTSEQ')[0]) > 0:
                    inlstseq = station['LSTSEQ'][()]
                else:
                    inlstseq = station['lstseq'][()]
                minlstseq[jjj] = np.min(inlstseq)
                maxlstseq[jjj] = np.max(inlstseq)
            except KeyError:
                f.close()
        except IOError:
            {}
    
    idx = np.searchsorted(maxlstseq, lstseq) 
    
    #empty arrays to be filled later
    x_target = np.zeros(len(idx))
    y_target = np.zeros(len(idx))
    error_mask=np.zeros(len(idx),dtype="bool")
    
    uidx = np.unique(idx)
    counter = 0
    
    for iii in uidx:
        curset = np.where(idx == iii)[0]
        fast = h5py.File(fastfiles[iii], "r")
        
        for jjj in curset:
            xx,yy, err = conv(RA_target, DEC_target,
                              lstseq[jjj],
                              lst[jjj],
                              jd[jjj],
                              fastfiles[iii], fast, camstr[camera])
            error_mask[counter]=err
            
            if err:
                x_target[counter], y_target[counter] = xx, yy
                with open("{}/x_y_{}_{}.txt".format(dir,star,camera), "a+") as file:
                    file.write("{}\t{}\t{}\n".format(lstseq[jjj],x_target[jjj],y_target[jjj]))
            else:
                x_target[counter],y_target[counter]=-1,-1
            counter = counter+1
        fast.close()
        
    #if any conversion did not succeed this is discarded from the data
    x_target=np.array(x_target)[error_mask]
    y_target=np.array(y_target)[error_mask]
    lstseq=lstseq[error_mask]
    lst=lst[error_mask]
    jd=jd[error_mask]
    files=files[error_mask]
    xref=xref[error_mask]
    yref=yref[error_mask]
    
    return (x_target, y_target,files,lstseq,lst,jd,xref,yref)
        
    
        
#make a post stamp image around the target        
def cube(filename,jd,xsel,ysel,R,plot,dir,x_ref,y_ref,star,t,camera,ASCC):
    """cube() takes a region around the target from an image and returns this.
INPUT:
filename: name of the fits file containing the original image
jd: the Julian date of the image
xsel: the x coordinate of the target of the image
ysel: the y coordinate of the target of the image    
R: the values determining the dimensions of the post stamp images (2Rx2R)
plot: Boolean value determining wheter a plot of the post stamp image is made
dir: string containing the path to the working directory
x_ref: the x coordinates of the reference stars of the image
y_ref: the y coordinates of the reference stars of the image
star: name of the target star
t: index of the files from the total array of files
camera: index of the camera as per burp.cams
ASCC: array with ASCC names of the reference stars
OUTPUT:
selection: the post stamp image around the target of the image    
    """
    
    image=fits.getdata(filename)
    xmin=int(round(xsel-R))
    xmax=int(round(xsel+R))
    ymin=int(round(ysel-R))
    ymax=int(round(ysel+R))
    
    selection=image[ymin:ymax,xmin:xmax]
    
    #if the selected region is too close to the edge to not be fully displayed
    #return region with zeros
    if xsel<R:
        print ("pech")
        return (np.zeros((2*R,2*R)))
    elif xsel+R>image.shape[1]:
        print ("pech")
        return (np.zeros((2*R,2*R)))
    if ysel<R:
        print ("pech")
        return (np.zeros((2*R,2*R)))
    elif ysel+R>image.shape[0]:
        print ("pech")
        return (np.zeros((2*R,2*R)))
    
    camstr = ['AUE', 'AUW', 'LSC', 'LSE', 'LSN', 'LSS', 'LSW', 'SAE', 'SAW', 'LPC', 'LPE', 'LPN', 'LPS', 'LPW']
    
    #plot the post stamp image        
    if plot==True:
        start_date=datetime.datetime(1858,11,17,0,0,0)
        date=start_date+datetime.timedelta(jd-2400000.5)
        date=date.strftime("%Y %b %d %H:%M:%S")
        
        plt.figure()
        im=plt.imshow(selection,norm=LogNorm())#,vmax=10000,vmin=1000)
        x=np.append(xsel,x_ref)
        y=np.append(ysel,y_ref)
        colours=["red","cyan","orange","lawngreen"]
        labels=[star.replace("_"," "),"ASCC {}".format(ASCC[0]),"ASCC {}".format(ASCC[1]),"ASCC {}".format(ASCC[2])]
        circles=np.zeros(len(x)).tolist()
        for i in range(len(x)):
            aper=plt.Circle((x[i]-xsel+R,y[i]-ysel+R),4.5,color=colours[i],fill=False)
            plt.gca().add_artist(aper)
            sky1=plt.Circle((x[i]-xsel+R,y[i]-ysel+R),6,color=colours[i],fill=False,ls="--")
            plt.gca().add_artist(sky1)
            sky2=plt.Circle((x[i]-xsel+R,y[i]-ysel+R),21,color=colours[i],fill=False,ls="--")
            plt.gca().add_artist(sky2)
            circles[i]=aper
        plt.title("{} with camera {} on {} UT".format(star.replace("_"," "), camstr[camera],date))
        
        plt.legend(circles,labels)
        plt.axis('off')
        cbar=plt.colorbar(im)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('ADU', rotation=270)
        plt.gca().invert_yaxis()
        plt.savefig("{}/stars_{}_{}.pdf".format(dir,star,t))
        plt.close()
        
    return (selection)
        
#detrend moet er nog uit, maar eerst testen of ie zo werkt
#getting the flux for the target star and reference stars
def find_flux(cube,camera,X,Y,star,dir):
    """find_flux() uses the code Remko provided photometry.py and 
configuration.py to find the flux of the target and the reference stars in the 
images.
INPUT:
cube: 3D array with the first dimension the amount of images and the other x and
      y containing the post stamp images
camera: index of the camera as per burp.cams
X: the x coordinates of the desired stars from the original images
Y: the y coordinates of the desired stars from the original images
star: name of the target star
dir: string containing the path to the working directory
OUTPUT:
Flux: 2D array containing the fluxes of the stars (column) per image (row) 
      without the sky
eFlux: 2D array containing the error on fluxes of the stars (column) per image 
       (row) with the sky error as well
Sky: 2D array containing the fluxes of the sky around the stars (column) per 
     image (row)
eSky: 2D array containing the error on the fluxes of the stars (column) per 
      image (row)
Peak: 2D array containing the peak fluxes of the stars (column) per image (row)
Pflag: 2D array containing the flags of the stars (column) per image (row)

            The flags take the following values.
        0   : All went well.
        1   : The star was too close to the edge of the image no photometry was 
              performed.
        2   : There was a pixel in the sky annulus with a value greater than 
              badpixval.
        4   : The sky value was negative.
        8   : The peak value was greater than badpixval.
        16  : One of the flux measurements was negative.
    """

    #empty arrays to be filled later
    Flux=[]
    eFlux=[]
    Sky=[]
    eSky=[]
    Peak=[]
    Pflag=[]
    
    #to clear the files
    open("{}/flux_{}_{}.txt".format(dir,star,camera),"w").close()
    open("{}/e_flux_{}_{}.txt".format(dir,star,camera),"w").close()
    open("{}/flag_flux_{}_{}.txt".format(dir,star,camera),"w").close()

    #determine the photometry properties
    phot = photometry.Photometry(cfg.aper_fast, cfg.skyrad)
    
    for i,el in enumerate(range(cube.shape[0])):
        #to compensate for the cutting of the cube
        xsel=X[i]-(X[i][0]-cube.shape[1]*0.5)
        ysel=Y[i]-(Y[i][0]-cube.shape[1]*0.5)
        
        #necessary fot phot.get_phot
        Xsel=np.array(xsel,dtype="float")
        Ysel=np.array(ysel,dtype="float")
        
        image=cube[i]
        flux, eflux, sky, esky, peak, pflag = phot.get_phot(image, Xsel, Ysel)
        
        with open("{}/flux_{}_{}.txt".format(dir,star,camera), "a+") as file:
            file.write("{}\t{}\t{}\t{}\n".format(flux[:,1][0],flux[:,1][1],flux[:,1][2],flux[:,1][3]))
        with open("{}/e_flux_{}_{}.txt".format(dir,star,camera), "a+") as file:
            file.write("{}\t{}\t{}\t{}\n".format(eflux[:,1][0],eflux[:,1][1],eflux[:,1][2],eflux[:,1][3]))
        with open("{}/flag_flux_{}_{}.txt".format(dir,star,camera), "a+") as file:
            file.write("{}\t{}\t{}\t{}\n".format(pflag[0],pflag[1],pflag[2],pflag[3]))
            
        Flux.append(flux[:,1])
        eFlux.append(eflux[:,1])
        Sky.append(sky)
        eSky.append(esky)
        Peak.append(peak)
        Pflag.append(pflag)
        
    
    #create arrays containing the outputs with the columns the different stars
    Flux=np.reshape(Flux,(cube.shape[0],-1))
    eFlux=np.reshape(eFlux,(cube.shape[0],-1))
    Sky=np.reshape(Sky,(cube.shape[0],-1))
    eSky=np.reshape(eSky,(cube.shape[0],-1))
    Peak=np.reshape(Peak,(cube.shape[0],-1))
    Pflag=np.reshape(Pflag,(cube.shape[0],-1))
    

    return (Flux,eFlux,Sky,eSky,Peak,Pflag)


#to make a mask that excludes the period of the nova eruption
def make_mask(data,before,after):
    """make_mask() makes a mask that excludes the part of the array between
the values given.
INPUT:
data: array with the data
before: start value of the excluded region
after: end value of the excluded region
OUTPUT:
mask: a mask with the indices apart from the indices between before and after
    """
    mask=np.where(data<before)[0]
    mask=np.append(mask,np.where(data>after)[0])
    return (mask)

#to filter out the periodicity in the given period
def periodicy(jd,period,mag1,mag2,mag3,mask):
    """periodicity() filters out the periodicity of the given period for the 
magnitudes obtained with all three reference stars. It uses the mask made with
make_mask() to exclude the period around the nova eruption.
INPUT:
jd: array containing the Julian dates
period: the period over which the periodicities are filtered
mag1: array of the magnitude of the target star obtained with reference star 1
mag2: array of the magnitude of the target star obtained with reference star 2
mag3: array of the magnitude of the target star obtained with reference star 3
mask: mask containing the indices of the data points excluding the nova eruption
      period
OUTPUT:
mag1_new: array of the magnitude of the target star obtained with reference star 1
          with the periodicity over {period} filtered out
mag2_new: array of the magnitude of the target star obtained with reference star 2
          with the periodicity over {period} filtered out
mag3_new: array of the magnitude of the target star obtained with reference star 3
          with the periodicity over {period} filtered out
    """

    #creating empty arrays to fill later
    mag1_new=np.zeros(len(jd))
    mag2_new=np.zeros(len(jd))
    mag3_new=np.zeros(len(jd))

    #the periods are divided into 100 equal size parts to find the periodicities
    for i in range(100):
        mask1=np.where((jd%period>=period*0.01*i)&(jd%period<period*0.01*(i+1)))[0]

        #apply the first mask to exclude the nova part in the calibration
        #but apply the found calibration to it anyway
        mask2=np.unique(np.append(mask,mask1))

        #sigma_clip() rejects outliers
        sa1=sigma_clip(mag1[mask2],sigma=1,iters=3) 
        sa2=sigma_clip(mag2[mask2],sigma=1,iters=3)
        sa3=sigma_clip(mag3[mask2],sigma=1,iters=3)

        mag1_new[mask1]=np.array(mag1[mask1])-np.mean(sa1)
        mag2_new[mask1]=np.array(mag2[mask1])-np.mean(sa2)
        mag3_new[mask1]=np.array(mag3[mask1])-np.mean(sa3)
        
    return (mag1_new,mag2_new,mag3_new)

#to apply the full calibration
def calibration(lstseq,jd,periods,mag1,mag2,mag3,mask):
    """calibration() applies the full calibration of all the desired peeriods
by applying periodicity(). NOTE: the first period is assumed to concern the
lstseq, whereas the others are assumed to concern the Julian date.
INPUT:
lstseq: array containing the lstseq
jd: array containing the Julian dates
periods: array containing the desired periods over which the periodicities are 
         filtered out in order of application
mag1: array of the magnitude of the target star obtained with reference star 1
mag2: array of the magnitude of the target star obtained with reference star 2
mag3: array of the magnitude of the target star obtained with reference star 3
mask: mask containing the indices of the data points excluding the nova eruption
OUTPUT:
mag1_new: 2D array containing the filtered out magnitude of the target star
          obtained with reference star 1 with the columns corresponding to the
          filtered out period that was indicated in {periods} cummulatively
mag2_new: 2D array containing the filtered out magnitude of the target star
          obtained with reference star 2 with the columns corresponding to the
          filtered out period that was indicated in {periods} cummulatively
mag3_new: 2D array containing the filtered out magnitude of the target star
          obtained with reference star 3 with the columns corresponding to the
          filtered out period that was indicated in {periods} cummulatively
    """
    #creating empty arrays to fill later
    mag1_new=np.zeros((len(mag1),len(periods)))
    mag2_new=np.zeros((len(mag2),len(periods)))
    mag3_new=np.zeros((len(mag3),len(periods)))
    
    #apply the filtering for each given period
    #each next filtering uses the filtered magnitude of the previous step
    #creating a cummulatively filtered magnitude
    for i in range(len(periods)):
        #the first period is assumed to concern the lstseq
        if i==0:
            mag1_new[:,i],mag2_new[:,i],mag3_new[:,i]=periodicy(lstseq,
                                            periods[i],mag1,mag2,mag3,mask)
        #the other periods are assumed to concern the Julian date
        else:
            mag1_new[:,i],mag2_new[:,i],mag3_new[:,i] = periodicy(jd,periods[i],
                           mag1_new[:,i-1],mag2_new[:,i-1],mag3_new[:,i-1],mask)
            
    return (mag1_new,mag2_new,mag3_new)


