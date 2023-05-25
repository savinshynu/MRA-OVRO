import numpy as np
import glob
import os

#path of image files for a single integration (multiple bands)
path = "/leo/savin/ovro-lwa/low_res/nfcor/mra4b/conv_4c*"


#Change here as well 
path_conv = "/leo/savin/ovro-lwa/low_res/nfcor/mra4b/conv_4c*done.image*"



#Getting image files
files = sorted(glob.glob(path))
print(files)

#opening the first or low frequency image to get the beam shape
ia.open(files[0])
beam = ia.restoringbeam()
print("Low frequency beam pars:")
print(beam)
ia.close()


#Convolve each high frequency image with a low frequency gaussian 2d beam

print("Convolving data")

for i,filename in enumerate(files):
    print(filename)
    outsplit = os.path.splitext(filename)[0]
    outname = outsplit+'_done.image'
    if i == 0:
       os.system("cp -r %s %s " % (filename, outname))
    else:           
        ia.open(filename)
        ia.convolve2d(outfile = outname, major = beam["major"], minor = beam["minor"], pa = beam["positionangle"], targetres = True)
        ia.close()

#Defining the box region, Change it for each MRA as they are different
#Add the box for each MRA 
box2 = rg.box([1600,1700],[2600,2700])
box5 = rg.box([1700,1800],[2500,2600])
box1 = rg.box([1600,1525],[2400,2325])
box3 = rg.box([1600,1850],[2400,2650])
box4 = rg.box([1700,1700],[2300,2300])


#Specify the box here
box = box4

files_conv = sorted(glob.glob(path_conv))

for filename in files_conv:
    print(filename)
    outsplit = os.path.splitext(filename)[0]
    outname = outsplit+'_sub.image'
    print("Final beam resolution")

    ia.open(filename)
    print(ia.restoringbeam())

    #Specifiy the box here
    ia.subimage(outfile = outname, region = box)
    ia.close()


