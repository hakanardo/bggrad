import os
from cffi import FFI
ffi = FFI()

mydir = os.path.join(os.path.dirname(__file__), 'maxflow-v3.02.src')
ffi.cdef(open(os.path.join(mydir, 'capi.h')).read())

lib = ffi.set_source("bggrad._maxflow", '#include "capi.h"',
                 include_dirs=[mydir],
                 extra_compile_args=['-Ofast', '-mtune=native', '-march=native'],
                 sources=['maxflow-v3.02.src/graph.cpp', 'maxflow-v3.02.src/maxflow.cpp', 'maxflow-v3.02.src/capi.cpp'],
                 )

if __name__ == "__main__":
    ffi.compile()
