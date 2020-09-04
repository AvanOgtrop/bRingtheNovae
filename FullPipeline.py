################################################################################
#Code developed by Aniek van Ogtrop                                            #
#Last edit: 2020 August                                                        #
#The used functions can be found in functionsFullPipeline.py                   #
################################################################################

#necessary imports
import os
from astropy.io import fits
from astropy.table import Table
import numpy as np
import time

import sys
sys.path.append('/disks/web1/bring/products/')
sys.path.append('/net/belterwijde/data2/stuik/code')
sys.path.append("/net/belterwijde/data2/vogtrop")
import functionsFullPipeline as ffp

import burp

"""
This file contains the full pipeline to get the light curve of a star when only
the RA and DEC are known.

0. Before the code can run some parameters need to be specified by the user,
   namely:
   RA_target (in degrees) (can use function ffp.get_to_degrees() to convert)
   DEC_target (in degrees) (can use function ffp.get_to_degrees() to convert)
   star (name of the target)
   camera (the camera used)
   directory (the working directory where everything is saved to)
   keep_text_files (boolean indicating if the txt files should be saved or removed)
   
1. First, the three closest stars that are not variable and have a blend of less
   than 0.1 are found as reference stars. 
   
2. Then the images that contain all three reference stars are found using the
   burp.detrend() method. 

3. Looping over these figures the RA and DEC of the target star is converted to
   x and y coordinates. 

3a.The failed conversions need to be filtered out. 

4. A data cube is made of the images to make the manipulating quicker later 
   on. 

5. With the x and y the flux can be found. 

6. The magnitude is found.

7. The data is saved in a fits table.
"""

#the catalogue that's used for the identifying of the stars
starcat=fits.open('/net/belterwijde/data2/bring/Catalogues/bringcat20180428.fits')[1].data

#only allowing stars with a magnitude of 8.4 or lower to be in the set
#as only lightcurves are available for those
starcatsel=starcat[~np.isnan(starcat["Vmag"])]
starcatsel=starcatsel[starcatsel["Vmag"]<=8.4]




#All the things that need to be specified to run the code
################################################################################
RA_target,DEC_target=ffp.get_to_degrees(10,36,15.43,-59,35,53.7) #degrees
star="Nova_Carinae_2018_TEST"
camera=0  #AUE=0, AUW=1, SAE=7, SAW=8

#get a new directory for every new star
basepath="/net/belterwijde/data2" 
user="vogtrop/Goal2" 
directory="{}_{}".format(star,camera)

#the data is saved in a fits Table, however in intermediate steps the parameters
#are also saved in separate txt files
#choose wheter to keep these or remove them
keep_text_files=False

#choose wheter to apply the calibration
calibration=True

if calibration==True:
    #Julian date marking the start of the nova eruption period to exclude this from
    #the filtering
    before=2458175
    #Julian date marking the end of the nova eruption period to exclude this from
    #the filtering
    after=2458290
    
    #the periods over which periodic trends are removed from the magnitudes in this
    #order
    periods=[13500,burp.sidereal,burp.synodic,1] 

################################################################################

#make a new directory (if it did not exist yet) to store the results from the 
#pipeline
dir = os.path.join(basepath,user,directory)
if not os.path.exists(dir):
    os.mkdir(dir)








    
################################################################################
#1. Find the three closest stars that meet the criteria                        #
################################################################################
    
idx,ASCC,HD,ra,dec,Vmag=ffp.find_3_closest(RA_target,DEC_target,starcatsel)

print ("\nreference stars found: {}, {}, {}.\n".format(ASCC[0],ASCC[1],ASCC[2]))
np.savetxt("{}/Vmag_{}_{}.txt".format(dir,star,camera),Vmag)
np.savetxt("{}/ASCC_{}_{}.txt".format(dir,star,camera),ASCC,fmt="%s")

################################################################################
#2. Get the detrend data from the 3 closest stars                              #
################################################################################

detrendStars=np.array([])
for i in range(len(idx)):
    detrend=burp.detrend(ASCC[i])
    detrendStars=np.append(detrendStars,detrend)

print ("\nthe detrend data is collected\n")

#get the names of the files and the x, y, jd and lstseq of these images of the
#reference stars
lstseq,x_ref,y_ref,jd,lst,files=ffp.get_lstseq(dir,detrendStars, camera, star)

print ("\nthe x and y of the reference stars are found and the files are saved\n")

################################################################################
#3. RA DEC conversion                                                          #
################################################################################

#get the x and y coordinates corresponding to the images
x_target,y_target,files,lstseq,lst,jd,x_ref,y_ref=ffp.RA_DEC_conversion(dir,RA_target,DEC_target,files,star,camera,lstseq,lst,jd,x_ref,y_ref)



print ("\nthe x and y of the target are found\n")

################################################################################
#4. Make a data cube                                                           #
################################################################################


#determining the R that should be used:
#the maximum distance between the target star and its reference stars is used.
#25 pixels are added to make sure the sky aperture of 21 pixels is included
max_dist_squared=np.max((np.array(x_ref)-np.array(x_target)[:,None])**2+(np.array(y_ref)-np.array(y_target)[:,None])**2)
R=np.sqrt(max_dist_squared)+25


#initialising the data cube
cubearray=np.zeros((len(files),2*int(round(R)),2*int(round(R))))

j=0
error_index=np.array([],dtype="int")
succeed_index=np.array([],dtype="int")

for i,el in enumerate(files):
    plot=False
    start=time.time()
    
    #make an image of the field every 500 files
    if i%500==0:
        plot=True       
        print (i,time.time()-start) 
        
    
    selection=ffp.cube(el,jd[i],x_target[i],y_target[i],int(round(R)),plot,dir,x_ref[i],y_ref[i],star,i,camera,ASCC)
    #when not the whole region is available because the region is too close to 
    #the edge of the image it is discarded
    if np.sum(selection)!=0:
        cubearray[j,:,:]=selection
        j+=1
        succeed_index=np.append(succeed_index,i)
    else:
        error_index=np.append(error_index,i)
    print (i,time.time()-start)

#if not all regions could be taken this needs to be corrected for in the other
#variables as well
if len(error_index)>0:
    cubearray=cubearray[0:len(succeed_index)]
    files=files[succeed_index]
    x_ref=x_ref[succeed_index]
    y_ref=y_ref[succeed_index]
    jd=jd[succeed_index]

    x_target=x_target[succeed_index]
    y_target=y_target[succeed_index]
    lstseq=lstseq[succeed_index]
    
    np.savetxt("{}/files_{}_{}.txt".format(dir,star,camera),files,fmt="%s")
    np.savetxt("{}/x_reference_{}_{}.txt".format(dir,star,camera),x_ref,fmt="%f")
    np.savetxt("{}/y_reference_{}_{}.txt".format(dir,star,camera),y_ref,fmt="%f")
    np.savetxt("{}/jd_{}_{}.txt".format(dir,star,camera),jd,fmt="%f")
    np.savetxt("{}/x_y_{}_{}.txt".format(dir,star,camera),np.array([lstseq,x_target,y_target]).T,fmt="%f")
    np.savetxt("{}/lstseq_{}_{}.txt".format(dir,star,camera),lstseq,fmt="%i")

#making the header for the data cube    
header=fits.Header()
header["CAMERA"] = burp.cams[camera],"Camera used for the observations"
header["STAR"] = star, "Name of the observed target star"
header["FILE"] = "{}/data_{}_{}.txt".format(dir,star,camera), "Name of the file containing all variables"
header["FILEFILE"] = "{}/files_{}_{}.txt".format(dir,star,camera), "Name of the file containing the file names of the original images"
fits.writeto('{}/Cube_{}_{}.fits'.format(dir,star,camera), cubearray,header, overwrite=True)

print ("\nthe data cube is made and saved\n")

############################################################################################
#5.+6. Obtaining the flux and magnitude                                                    #
############################################################################################

#first the x arrays and y arrays need to be combined
X=np.hstack((x_target.reshape((len(x_target),-1)),x_ref))
Y=np.hstack((y_target.reshape((len(y_target),-1)),y_ref))

flux, eflux, sky, esky, peak, flag=ffp.find_flux(cubearray,camera,X,Y,star,dir)


#converting the flux to magnitude
def mag(m2,F1,F2):
    """mag() converts the flux of the target to the magnitude using the flux and
magnitude of a reference star.
INPUT:
m2: magnitude of the reference star
F1: flux of the target
F2: flux of the reference star
OUTPUT:
magnitude of the target
    """
    return (m2-2.5*np.log10(F1/F2))

#calculating the uncertainty on the magnitude
def emag(eF1,eF2,F1,F2):
    """emag() finds the uncertainty in the magnitude of the target based on the
uncertainty of the flux of both the target and the reference star.
INPUT:
eF1: uncertainty on the flux of the target
eF2: uncertainty on the flux of the reference star
F1: flux of the target
F2: flux of the reference star
OUTPUT:
uncertainty on the magnitude of the target
    """
    return(np.sqrt(2.5**2/(np.log(10))**2*((eF1/F1)**2+(eF2/F2)**2)))


M1=mag(Vmag[0],flux[:,0],flux[:,1])
M2=mag(Vmag[1],flux[:,0],flux[:,2])
M3=mag(Vmag[2],flux[:,0],flux[:,3])

eM1=emag(eflux[:,0],eflux[:,1],flux[:,0],flux[:,1])
eM2=emag(eflux[:,0],eflux[:,2],flux[:,0],flux[:,2])
eM3=emag(eflux[:,0],eflux[:,3],flux[:,0],flux[:,3])

np.savetxt("{}/magnitude_{}_{}.txt".format(dir,star,camera),np.array([M1,M2,M3]).T,fmt="%f")
np.savetxt("{}/e_magnitude_{}_{}.txt".format(dir,star,camera),np.array([eM1,eM2,eM3]).T,fmt="%f")

np.savetxt("{}/sky_{}_{}.txt".format(dir,star,camera),sky,fmt="%f")
np.savetxt("{}/e_sky_{}_{}.txt".format(dir,star,camera),esky,fmt="%f")



print ("\nthe flux and magnitue are found\n")


######################################################################################
#7. Making a fits table of all information                                           #
######################################################################################

#making the labels of the columns
keys=["lstseq","jd","flux_target","flux_r1","flux_r2","flux_r3","e_flux_target",
    "e_flux_r1","e_flux_r2","e_flux_r3","sky_target","sky_r1","sky_r2","sky_r3",
    "e_sky_target","e_sky_r1","e_sky_r2","e_sky_r3","flag_target","flag_r1",
    "flag_r2","flag_r3","x_target","x_r1","x_r2","x_r3","y_target","y_r1",
    "y_r2","y_r3","M_r1","M_r2","M_r3","e_M_r1","e_M_r2","e_M_r3"]

#combining the different arrays    
values=[np.array(lstseq,dtype="int32"),np.array(jd,dtype="float64"),
 np.array(flux[:,0],dtype="float64"),np.array(flux[:,1],dtype="float64"),
 np.array(flux[:,2],dtype="float64"),np.array(flux[:,3],dtype="float64"),
 np.array(eflux[:,0],dtype="float64"),np.array(eflux[:,1],dtype="float64"),
 np.array(eflux[:,2],dtype="float64"),np.array(eflux[:,3],dtype="float64"),
 np.array(sky[:,0],dtype="float64"),np.array(sky[:,1],dtype="float64"),
 np.array(sky[:,2],dtype="float64"),np.array(sky[:,3],dtype="float64"),
 np.array(esky[:,0],dtype="float64"),np.array(esky[:,1],dtype="float64"),
 np.array(esky[:,2],dtype="float64"),np.array(esky[:,3],dtype="float64"),
 np.array(flag[:,0],dtype="int32"),np.array(flag[:,1],dtype="int32")
 ,np.array(flag[:,2],dtype="int32"),np.array(flag[:,3],dtype="int32"),
 np.array(X[:,0],dtype="float64"),np.array(X[:,1],dtype="float64"),
 np.array(X[:,2],dtype="float64"),np.array(X[:,3],dtype="float64"),
 np.array(Y[:,0],dtype="float64"),np.array(Y[:,1],dtype="float64"),
 np.array(Y[:,2],dtype="float64"),np.array(Y[:,3],dtype="float64"),
 np.array(M1,dtype="float64"),np.array(M2,dtype="float64"),
 np.array(M3,dtype="float64"),np.array(eM1,dtype="float64"),
 np.array(eM2,dtype="float64"),np.array(eM3,dtype="float64")]
    
#making the header for the table
header={"star": "{}".format(star),
        "camera": "{}".format(camera),
        "RA_targ": "{}".format(RA_target),
        "DEC_targ": "{}".format(DEC_target),
        "r1_name": "{}".format(ASCC[0]),
        "r2_name": "{}".format(ASCC[1]),
        "r3_name": "{}".format(ASCC[2]),
        "r1_Vmag": "{}".format(Vmag[0]),
        "r2_Vmag": "{}".format(Vmag[1]),
        "r3_Vmag": "{}".format(Vmag[2])}


#saving the table
table1=Table(values,names=keys,meta=header)
table1.write("{0}/data_{1}_{2}.fits".format(dir,star,camera),overwrite=True)


#remove the text files if preferred
if keep_text_files==False:
    os.remove("{}/ASCC_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/e_flux_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/e_magnitude_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/e_sky_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/files_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/flag_flux_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/flux_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/jd_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/lstseq_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/magnitude_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/sky_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/Vmag_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/x_reference_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/y_reference_{}_{}.txt".format(dir,star,camera))
    os.remove("{}/x_y_{}_{}.txt".format(dir,star,camera))

######################################################################################
#8. Filtering out periodicities                                                      #
######################################################################################

#create a mask that excluded the data point containing the nova eruption
mask=ffp.make_mask(jd,before,after)

#apply the filtering for the given periods
M1_new,M2_new,M3_new=ffp.calibration(lstseq,jd,periods,M1,M2,M3,mask)

#making the labels of the columns
keys=["lstseq","jd","M_r1_fil1","M_r1_fil12","M_r1_fil123","M_r1_fil1234","M_r2_fil1",
      "M_r2_fil12","M_r2_fil123","M_r2_fil1234","M_r3_fil1","M_r3_fil12","M_r3_fil123",
      "M_r3_fil1234"]

#combining the different arrays 
values=[np.array(lstseq,dtype="int32"),np.array(jd,dtype="float64"),
 np.array(M1_new[:,0],dtype="float64"),np.array(M1_new[:,1],dtype="float64"),
 np.array(M1_new[:,2],dtype="float64"),np.array(M1_new[:,3],dtype="float64"),
 np.array(M2_new[:,0],dtype="float64"),np.array(M2_new[:,1],dtype="float64"),
 np.array(M2_new[:,2],dtype="float64"),np.array(M2_new[:,3],dtype="float64"),
 np.array(M3_new[:,0],dtype="float64"),np.array(M3_new[:,1],dtype="float64"),
 np.array(M3_new[:,2],dtype="float64"),np.array(M3_new[:,3],dtype="float64")]

#making the header for the table
header={"star": "{}".format(star),
        "camera": "{}".format(camera),
        "RA_targ": "{}".format(RA_target),
        "DEC_targ": "{}".format(DEC_target),
        "r1_name": "{}".format(ASCC[0]),
        "r2_name": "{}".format(ASCC[1]),
        "r3_name": "{}".format(ASCC[2]),
        "r1_Vmag": "{}".format(Vmag[0]),
        "r2_Vmag": "{}".format(Vmag[1]),
        "r3_Vmag": "{}".format(Vmag[2])}

#saving the table
table1=Table(values,names=keys,meta=header)
table1.write("{0}/data_filtered_{1}_{2}.fits".format(dir,star,camera),overwrite=True)


