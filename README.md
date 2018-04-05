# generate_lightcurve

Readme coming soon.

Important to do:

Write nice documentation

Implement error bars with the calculation of intensity

Portability and quality of life to do:

automatically create a shapemodels directory if not already exist

for each asteroid, automatically create a shapemodels/asteroidnumber director if not already exist

Add a way for the code to  serch for and download all necessary shape and spin files from DAMIT

e.g. replacing the "3200" in

http://astro.troja.mff.cuni.cz/projects/asteroids3D/web.php?page=db_listing&mode=search&asteroid_id=&model_id=&asteroid_number=3200&asteroid_name=&asteroid_designation=

and automatically downloading:

http://astro.troja.mff.cuni.cz/projects/asteroids3D/data/archive/1-1000/A278.M399.shape.txt

http://astro.troja.mff.cuni.cz/projects/asteroids3D/php/spin.txt.php?model_id=399

have the code automatically fine location of oorb executable

mdfind -name oorb | grep oorb-master/main/oorb

have code find MPCOORB.DAT

mdfind -name "MPCORB.DAT"

have code find lcgenerator

/Users/bolin/NEO/lightcurves/version_0.2.1/lcgenerator/lcgenerator

#have code find working directory for 
the sys.path.insert(0, '/Users/bolin/NEO/lightcurves/scripts/')
from lightcurve_functions import *

pwd
#add list of requirements

linux:

mdfind
pwd

lcgenerator
http://astro.troja.mff.cuni.cz/projects/asteroids3D/web.php?page=download_software

oorb
https://github.com/oorb/oorb

download:

MPCORB.DAT

python:

string
os
numpy as np
re
argparse
sys
astropy
random