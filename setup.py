from setuptools import setup
from bggrad import __version__

setup(
    name='bggrad',
    description='Background/Foreground segmentation based on gradient directions',
    long_description='''
    Implementation of "Bayesian Formulation of Gradient Orientation Matching",
    https://link.springer.com/chapter/10.1007/978-3-319-20904-3_9
    ''',
    version=__version__,
    packages=['bggrad'],
    zip_safe=False,
    url='https://github.com/hakanardo/bggrad',
    author='Hakan Ardo',
    author_email='hakan@debian.org',
    license='MIT',
    setup_requires=["cffi>=1.0.0"],
    cffi_modules=["build_bggrad.py:ffi"],
    install_requires=["cffi>=1.0.0", "numpy>=1.7.1"],
)
