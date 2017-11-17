import os
from cffi import FFI
ffi = FFI()

mydir = os.path.abspath(os.path.dirname(__file__))
ffi.cdef(open("src/bggrad.h").read())
ffi.set_source("bggrad._bggrad", '#include "src/bggrad.h"',
               extra_compile_args=['-Ofast', '-mtune=native', '-march=native'],
               sources=["src/bggrad.c"],
               include_dirs=[mydir])

if __name__ == "__main__":
    ffi.compile()
