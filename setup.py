'''
from distutils.core import setup
import py2exe
import sys

setup(console=['main.py'])
'''
#setup.py


from cx_Freeze import setup, Executable

setup(
    name = "ideogram",
    version = "1.0.0",
    options = {"build_exe": {
        'packages': ["os","matplotlib"],
        'include_files': ['icon.png'],
        'excludes': ['collections.abc'],
        'include_msvcr': True,
    }},
    executables = [Executable("main.py",base=None)]
    )
