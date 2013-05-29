"""
NeuroPy
=======

analysis and modeling of electrophysiological data
for questions/comments contact <kilian@koepsell.de>
"""


#
# define package variables
#
import os
datapath = os.path.join(os.sep,'data')
cachepath = os.path.join(os.sep,'data','cache')


# We wrap all graphics-related stuff inside pybrain so we can run in
# environments with no access to any plotting (adapted from fernando perez)
try:
    import pylab
    has_matplotlib = True
except:
    has_matplotlib = False
    class Dummy:
        """Dummy do-nothing class.

        Instances of this class return a dummy attribute on all accesses, which
        can be called."""

        def __init__(self,name):
            print "*** WARNING *** Dummy object replacing %r." % name

        def __str__(self):
            return "<DUMMY-OBJECT>"

        __repr__ = __str__

        def __getattr__(self,name):
            return self.dummy

        def dummy(self,*args,**kw):
            """Dummy function, which doesn't do anything (returns None)."""
            return None

    pylab = Dummy('pylab')

# we wrap pystream so that we can use a GPU whenever present
try:
    import pystream
    print '*** found pystream ***'
    del pystream
    from pystream import cudart
    cudart.setDevice(0)
    print '*** found GPU ***'
    del cudart
    has_pystream = True
except:
    has_pystream = False

#
# system and user specific settings
#
import sys
__dir__ = os.path.dirname(__file__)
try:
    __user__ = os.getlogin()
except:
    __user__ = 'nobody'
try:
    __uid__ = str(os.getuid())
except:
    __uid__ = 'none'
__sys__ = sys.platform
__init_user = os.path.join(__dir__,'__init__'+__user__+'.py')
__init_system = os.path.join(__dir__,'__init__'+__sys__+'.py')
__init_system_user = os.path.join(__dir__,'__init__'+__sys__+'_'+__user__+'.py')
__init_system_uid = os.path.join(__dir__,'__init__'+__sys__+'_'+__uid__+'.py')

try:
    # print >> sys.stdout, 'loading init file '+__init_system
    execfile (__init_system)
except:
    # print >> sys.stdout, 'no init file '+__init_system
    pass

try:
    # print >> sys.stdout, 'loading init file '+__init_system_user
    execfile (__init_system_user)
except:
    # print >> sys.stdout, 'no init file '+__init_system_user
    pass

try:
    # print >> sys.stdout, 'loading init file '+__init_system_uid
    execfile (__init_system_uid)
except:
    # print >> sys.stdout, 'no init file '+__init_system_uid
    pass

#user should be loaded last, to overwrite whatever setting were set by the other files
try:
    print >> sys.stdout, 'loading init file '+__init_user
    execfile (__init_user)
except:
    print >> sys.stdout, 'no init file '+__init_user
    pass

#
# load modules
#
modules = ['utils','datalib','data','dsplib','spikelib','model','plotlib','hdf5lib', 'bdflib', 'eyelib']
# add custom modules
from fnmatch import fnmatch
modules += [m[:-3] for m in os.listdir(__dir__) if fnmatch(m,'custom_*.py')]
__all__ = modules[:]
gvar,lvar = globals(),locals()
for name in modules:
    mod = __import__(name,gvar,lvar,[])
    # reload modules (useful during development)
    reload(mod)
    __all__ += mod.__all__

for name in modules:
    exec('from %s import *'%name)

def test():
    from numpy.testing import NumpyTest
    NumpyTest().test()

# Namespace cleanup
del name,gvar,lvar,os,sys,fnmatch

# def test(verbosity=1):
#     """Run the entire mwadap test suite, at the given verbosity.
# 
#     Valid values of verbosity are 0,1,2."""
#
#     from test import test
#     test.run(verbosity)
