# Sparse Noise
from Dimstim.SparseNoise import *

# Static parameters, remain constant over entire experiment

# pre-experiment duration to display blank screen (sec)
preexpSec  = 1
# post-experiment duration to display blank screen (sec)
postexpSec = 1
# target ori offset (deg), also the grid ori, read from manbar 1 if left as None
orioff    = None

# grid width (number of cells)
ncellswide = 9
# grid height (number of cells)
ncellshigh = 9
# grid width (deg), is set to multiple of widthDeg[0] (no overlap) if left as None
regionwidthDeg  = 8*ncellswide/3 #1*ncellswide
# grid height (deg), is set to multiple of heightDeg[0] (no overlap) if left as None
regionheightDeg = 8*ncellswide/3 #None #0.333*ncellshigh

# Dynamic variables, could conceivably change between sweeps. If a variable is assigned multiple values in a [list], then its name must be entered as a variable in the varlist if you want it to vary across sweeps.
# Use span(start,stop,step,n) f'n to define start and 2 of: stop; step; number of values
# Use len(var) to get the # of values in a variable, use list*n to make n concatenated copies of list, handy for generating variables with the same # of values as some other existing variable.
# Use scale(list,inv=True,factor) to get the element-wise inverse of a list, and then optionally scale each element by a factor. Useful for making two variables inversely covary. Don't forget to make them part of the same dimension.

# index of cell horizontal position (origin at left)
xi = range(0,ncellswide)
# index of cell vertical position (origin at bottom)
yi = range(0,ncellshigh)
# target ori relative to grid ori == orioff (deg)
ori       = 0#[-15,0,15]#span(start=0,stop=360-18,step=18)
# target width (deg), read from manbar if left as None
widthDeg  = 0.333#[0.5,0.75]
widthDeg  = 0.210#[0.5,0.75]
widthDeg  = 0.090#[0.5,0.75]
# target height (deg), read from manbar if left as None
heightDeg = 0.333#[0.5,0.75]
heightDeg = 0.210#[0.5,0.75]
heightDeg = 0.090#[0.5,0.75]
# target brightness (0-1)
brightness   = [1]#[0,0.2,0.4,0.6,0.8,1]
brightness   = 1#[0,0.2,0.4,0.6,0.8,1]
# antialiase the target? (0=no, 1=yes), eg assign it [0,1,1]*n and/or shuffle it
antialiase   = 1
# background brightness (0-1)
bgbrightness = 0.5
bgbrightness   = [0.2,0.8]
bgbrightness   = [0,0.2,0.4,0.6,0.8,1]
bgbrightness = [0,0]*10
bgbrightness.extend(numpy.linspace(0.,1., 256)[:200].tolist())#span(start=0.0,stop=1.0,n=11)#span(start=0,stop=1,n=len(ori))
bgbrightness.extend(numpy.linspace(0.,1., 256)[200::-1].tolist())#span(start=0.0,stop=1.0,n=11)#span(start=0,stop=1,n=len(ori))
# antialiase the target? (0=no, 1=yes), eg assign it [0,1,1]*n and/or shuffle it
# sweep duration (msec)
sweeptimeMsec = 1000/20 / 2.56
# post-sweep duration to display blank screen (msec)
postsweepMsec = 0

# Each entry in 'varlist' represents one variable. Each entry is itself made up of a 2-entry dictionary with 'shuffle' and 'dim' fields. The 'shuffle' field is a flag (0=leave ordered, 1=shuffle, 2=randomize). shuffle==1 shuffles the variable's dimension during the experiment, shuffle==2 randomly samples it instead. Variable position in the varlist doesn't matter. However, the variable's associated 'dim' value relative to all the other 'dim' values in the variable list determines its order in the nested for loops that generate the combinations of values for each sweep: the variable with the lowest 'dim' value becomes the outermost for loop and changes least often; the variable with the highest 'dim' value becomes the innermost for loop and changes on every sweep. 'dim' must be an integer (+ve, -ve, or 0). Variables with the same 'dim' value are part of the same dimension, are shuffled together, and must therefore be assigned the same number of values and the same shuffle flag. varlist should probably never be left empty.

varlist = {
           'xi':           {'shuffle':0, 'dim':1},
           'yi':           {'shuffle':0, 'dim':2},
#           'ori':          {'shuffle':1, 'dim':1},
#           'widthDeg':     {'shuffle':1, 'dim':4},
#           'heightDeg':    {'shuffle':1, 'dim':5},
           #'brightness':   {'shuffle':1, 'dim':4},
#           'antialiase':   {'shuffle':0, 'dim':5},
           'bgbrightness' :{'shuffle':0, 'dim':3},
#           'sweeptimeMsec':{'shuffle':1, 'dim':5},
#           'postsweepMsec':{'shuffle':0, 'dim':3},
          }

# number of times to run the whole combination of variables
nruns              = 1
# shuffle/randomize variables across runs? (0=no, runs are identical, 1=yes, runs are different)
shuffleRuns        = 1
# (n,duration) do a blank sweep every n'th sweep for duration in seconds, n must be in: {0=off ,2,3,4....}
blankSweep         = (0,0)
# shuffle the position of the blank sweeps? (0=no, 1=yes)
shuffleBlankSweeps = 0
colour=[1,0,0]

###########################################################################################
# Run it
sparsenoise(preexpSec, postexpSec, orioff,
            ncellswide, ncellshigh, regionwidthDeg, regionheightDeg,
            xi, yi, ori, widthDeg, heightDeg, brightness, colour, antialiase, bgbrightness, sweeptimeMsec, postsweepMsec,
            varlist, nruns, shuffleRuns, blankSweep, shuffleBlankSweeps)
