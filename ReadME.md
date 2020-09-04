# bRing the Novae

This folder contains the code that has been used in the bRing the Novae project.

The code allows the user to extract light curves of any given right ascension and declination within the field of view of bRing. The light curve can be detrended for different periodicities indicated by the user.

The code is developed by Aniek van Ogtrop and is last edited in September 2020.

## Usage

Both python files are extensively commented so they should be straightforward to use. Both files should be downloaded. FullPipeline.py contains the working file and functionsFullPipeline.py contains all the functions associated with FullPipeline.py. Please note that the code is made to work in python 2 and does not work in python 3.

FullPipeline.py shows where the user should specify its preferences, but for clarity it is reiterated here as well:

The first thing that needs to be specified is the RA and DEC of the target. This needs to be in degrees, however a function is made to convert hour notation to degrees; ffp.get_to_degrees().

```python
RA_target,DEC_target=ffp.get_to_degrees(10,36,15.43,-59,35,53.7) #degrees
```

The second thing to specify is the name of the target which is used in plots and the names of files.

```python
star= "target"
```

The preferred camera of bRing also needs to be specified. The cameras are indicated with an index where AUE=0, AUW=1, SAE=7 and SAW=8.

```python
camera=0 #AUE=0, AUW=1, SAE=7, SAW=8
```

The directory where the files are saved needs to be specified.

```python
directory="/directory"
```

During intermediate steps results are saved in individual txt files. The user can choose to keep these or remove them.

```python
keep_text_files=True
```

The code can also calibrate out periodicities in the light curves. The user can choose to implement this.

```python
calibration=True
```

If the calibration is implemented the periods over which periodicities are found need to be specified. Please note that the code assumes that the first period concerns the lstseq and the next concern the Julian date.

```python
periods=[13500,burp.sidereal,burp.synodic,1] 
```

The periodicities are found in the period that does not contain the nova eruption. Hence the period of the nova eruption should be specified. This is done in Julian date.

```python
before=2458175 #indicating the start of the nova eruption period
after=2458290  #indicating the end of the nova eruption period
```

After the specifications of these parameters, the code will work as it should.