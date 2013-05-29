"""
This module contains utilities that don't fit any category
"""

import numpy as np
import sys,os,re,time

__all__ = ['flush','make_dir','now','tic','toc','save_cache','load_cache','round_up','smod','ipshell','ipython',
           'nchoosek','ksubset','ksubsets','eps']
eps = np.MachAr().eps

def flush():
    sys.stdout.flush()

def ipshell():
    '''
    starts an embedded IPython shell
    '''
    from IPython.Shell import IPShellEmbed
    args = ['-pylab', '-colors','NoColor']
    ipshell = IPShellEmbed(args, banner = 'Dropping into IPython...', exit_msg = 'Leaving IPython...')
    return ipshell()

def ipython():
    '''
    replaces current session by a full IPython session
    '''
    from IPython.Shell import start
    return start().mainloop()

def make_dir(pathname):
    if not os.path.exists(pathname):
        try:
            # makedirs is like mkdir -p, creates all subdirectories, not just the right most one (was os.mkdir)
            os.makedirs(pathname) 
        except:
            if not os.path.exists(pathname):
                assert 0, 'cannot create directory '+pathname


if sys.platform == 'win32':
    now = time.clock
else:
    now = time.time


_timer = now()
def tic(timer='timer'):
    print 'tic (%s) ... '%timer
    globals()['_'+timer] = now()
def toc(timer='timer'):
    _timer = globals().get('_'+timer,'_timer')
    dt = now()-_timer
    if dt < 1: print "toc (%s) ... %.1f ms"% (timer,1000.*dt)
    else: print "toc (%s) ... %.1f s"% (timer,dt)
    return dt

##
## caching of computed data
##
def load_cache(arrayname,cachedir,dtype='d'):
    """
    load pre-computed data from cache file
    """
    filename = os.path.join(cachedir,'%s.cache'%arrayname)
    try:
        return np.fromfile(filename,dtype=dtype)
    except:
        return None

def save_cache(data,arrayname,cachedir):
    """
    save computed array data to cache file
    """
    make_dir(cachedir)
    data.tofile(os.path.join(cachedir,'%s.cache'%arrayname))


def parse_dir(path,prefix=None,suffix=None,key_type=int):
    if suffix != None: # remove trailing dot
        suffix = re.sub('^\.','',suffix)
    if prefix == None:
        p = re.compile('(\w+)\.'+suffix+'$')
    else:
        p = re.compile('^'+prefix+'(\w+)\.'+suffix+'$')
    filelist = os.listdir(path) # list of all files in directory
    filenames = {}              # dict of matching files with id
    for fn in filelist:
        match = p.match(fn)
        if match == None: continue
        filenames[key_type(match.group(1))] = match.group(0)
    return filenames

def round_up(x):
    """
    Rounds the number x to the next number with one significant digit

    The return value has the form m*10**n with integers m,n
    """
    m = 10**np.floor(np.log10(abs(x))) or 1
    return m*(np.floor(x/m)+1)

def smod(x,p=2*np.pi):
    "Returns x mod p, symmetrically centered at zero"
    return np.mod(x+p/2,p)-p/2.


def nchoosek(n,k):
    """
    nchoosek(n,k) returns n choose k
    """
    if k > n:
        return 0
    if k == 0:
        return 1
    if n <= 1:
        return 1
    else:
        return nchoosek(n-1,k)+nchoosek(n-1,k-1)

class ksubset:
    """
    ksubset(n,k) returns a generator of all k-element subsets of range(n)

    see also ksubsets(n,k)
    """
    def __init__(self,n,k):
        self.count = nchoosek(n,k)
        self.n = n
        self.k = k
        self.reset()
    def reset(self):
        self.__index = 0
        self.__current = self.subset(self.n,self.k,0)
    def __iter__(self):
        return self
    def __getitem__(self,ind):
        return self.subset(self.n,self.k,ind)
    def next(self):
        if self.__index == self.count:
            self.reset()
            raise StopIteration
        else:
            self.__current = self.subset(self.n,self.k,self.__index)
            self.__index += 1
        return self.__current
    def subset(self,n,k,m):
        f=0
        if (k==0)|(m>=nchoosek(n,k)):
            return []
        def ntot(f):
            return nchoosek(n-f-1,k-1)
        while f<n-k:
            if m < ntot(f):
                subs = [f]
                subs.extend(np.array(self.subset(n-f-1,k-1,m))+f+1)
                return subs
            else:
                m = m-ntot(f)
                f = f+1
        return range(n-k,n)

def ksubsets(n,k):
    """
    ksubsets(n,k) returns all n-choose-k k-element subsets of range(n)

    example:
    >>> nchoosek(4,2)
    [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]

    """
    return [i for i in ksubset(n,k)]
