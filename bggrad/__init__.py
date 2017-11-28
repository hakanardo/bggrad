try:
    from bggrad.api import BgGrad
except ImportError:
    print("Use build_* scripts to build C extention or use setup.py to install.")

__version_info__ = (1, 0, 0)
__version__ = '.'.join(str(i) for i in __version_info__)
