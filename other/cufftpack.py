""" Testing a more Pythonic interface to the FFT libraries.

NOTE: works, but segfaults on exit.
"""

import atexit
from ctypes import POINTER, c_float

import numpy as np

from pystream.cudaarray import RawCudaArray, CudaArrayFromArray
from pystream import cudart, cufft


c_complex = c_float*2

class PlanCache(object):
    """ Simple object that maintains a certain number of plans.

    Old plans get destroyed.
    """

    def __init__(self, size=64):
        # XXX: is that a good size?
        self.size = size

        # The map of plan arguments to plan handles and back again.
        self.handle_map = {}
        self.handle_map_inverse = {}

        # The list of handles in order of most recent use.
        self.handles = []

        # Make sure we clean up before shutting down the program.
        atexit.register(self.cleanup)

    def lookup(self, dims, type_, batch=None):
        """ Look up a plan in the cache or create one.

        If we need to create a new plan and the cache is full, dump the one
        least-recently used.
        """
        plan = self.handle_map.get((dims, type_, batch), None)
        if plan is None:
            plan = self.createPlan(dims, type_, batch)
        else:
            # Freshen the list of handles to ensure we're on top.
            i = self.handles.index(plan)
            self.handles.insert(0, self.handles.pop(i))

        return plan

    def createPlan(self, dims, type_, batch=None):
        """ Create a new plan.
        """
        ndims = len(dims)
        if ndims not in (1, 2, 3):
            raise ValueError("only 1, 2, and 3 dimensions are supported")
        if batch is not None and ndims != 1:
            raise ValueError("batching is only supported with 1D FFTs")
        if ndims == 1:
            args = (dims[0], type_, batch)
        else:
            args = dims + (type_,)
        func = {
            1: cufft.cufftPlan1d,
            2: cufft.cufftPlan2d,
            3: cufft.cufftPlan3d,
        }[ndims]
        plan = func(*args)
        self.handle_map[args] = plan
        self.handle_map_inverse[id(plan)] = args
        self.handles.insert(0, plan)

        if len(self.handles) > self.size:
            # Pop a handle.
            old = self.handles.pop()
            oldargs = self.handle_map_inverse.pop(id(old))
            self.handle_map.pop(oldargs)
            cufft.cufftDestroy(old)

        return plan

    def cleanup(self):
        """ Destroy all of the handles in the cache.
        """
        self.handle_map.clear()
        self.handle_map_inverse.clear()
        while self.handles:
            plan = self.handles.pop()
            cufft.cufftDestroy(plan)


_plan_cache = PlanCache()


def fft(a, out=None):
    """ Do a 1D forward FFT.
    """
    a = np.ascontiguousarray(a)
    if np.issubdtype(a.dtype, np.complexfloating):
        a = a.astype(np.csingle)
        type_ = cufft.CUFFT_C2C
        ct = c_complex
    else:
        a = a.astype(np.single)
        type_ = cufft.CUFFT_R2C
        ct = c_float
    nx = a.shape[-1]
    if a.ndim == 1:
        batch = None
    elif a.ndim == 2:
        # Batch up many 1D transforms.
        batch = a.shape[0]
    else:
        raise ValueError("only 1and 2 dimensions are supported")
    plan = _plan_cache.lookup((nx,), type_, batch)

    if out is None:
        out = np.empty(a.shape, np.csingle)
    else:
        if not isinstance(out.dtype , np.csingle):
            raise ValueError("output must be single-precision complex")
        if out.shape != a.shape:
            raise ValueError("output must have the same shape as the input")

    # Since we're always copying to memory on the card, we can always do
    # C2C transforms in-place!
    input = CudaArrayFromArray(a)
    if type_ == cufft.CUFFT_R2C:
        output = RawCudaArray(a.size, np.csingle)
        args = (plan, input.data, output.data)
        func = cufft.cufftExecR2C
    else:
        output = input
        args = (plan, input.data, input.data, cufft.CUFFT_FORWARD)
        func = cufft.cufftExecC2C

    func(*args)
    cudart.threadSynchronize()

    output.toArray(out)
    input.free()
    output.free()

    return out

def ifft(a, out=None):
    """ Do a 1D inverse FFT.
    """
    a = np.ascontiguousarray(a)
    a = a.astype(np.csingle)
    type_ = cufft.CUFFT_C2C
    ct = c_complex
    nx = a.shape[-1]
    if a.ndim == 1:
        batch = None
    elif a.ndim == 2:
        # Batch up many 1D transforms.
        batch = a.shape[0]
    else:
        raise ValueError("only 1and 2 dimensions are supported")
    plan = _plan_cache.lookup((nx,), type_, batch)

    if out is None:
        out = np.empty(a.shape, np.csingle)
    else:
        if not isinstance(out.dtype , np.csingle):
            raise ValueError("output must be single-precision complex")
        if out.shape != a.shape:
            raise ValueError("output must have the same shape as the input")

    # Since we're always copying to memory on the card, we can always do
    # C2C transforms in-place!
    input = CudaArrayFromArray(a)
    args = (plan, input.data, input.data)
    func = cufft.cufftExecC2C
    args += (cufft.CUFFT_INVERSE,)

    func(*args)
    cudart.threadSynchronize()

    input.toArray(out)
    input.free()
    input.free()

    return out


if __name__ == '__main__':
    import sys
    import time

    if sys.platform == 'win32':
        now = time.clock
    else:
        now = time.time

    cudart.setDevice(0)
    n = 2 ** 18
    batch = 16
    a = np.random.random((batch,n)).astype(np.csingle)

    t0 = now()
    f = fft(a)
    a2 = ifft(f) / n  # must normalize ourselves
    dt = now() - t0
    print 'CUDA FFT completed in %s s.' % dt
    print 'L1 error: %s' % np.absolute(a2 - a).sum()

    t0 = now()
    f = np.fft.fft(a)
    a2 = np.fft.ifft(f)
    dt = now() - t0
    print 'numpy FFT completed in %s s.' % dt
    print 'L1 error: %s' % np.absolute(a2 - a).sum()

