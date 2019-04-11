'''
from distutils.core import setup
import py2exe
#import sys

setup(console=['main.py'])

#setup.py

'''

'''
from cx_Freeze import setup, Executable

setup(
    name = "Astyanax,",
    version = "1.0.0",
    options = {"build_exe": {
        'packages': ["os","matplotlib"],
        'include_files': ['icon.png'],
        'excludes': ['collections.abc'],
        'include_msvcr': True,
    }},
    executables = [Executable("main.py",base=None)]
    )

'''

'''
import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["os", 'matplotlib'], "excludes": ["tkinter"]}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "Astyanax",
        version = "0.1",
        description = "My GUI application ideogram!",
        options = {"build_exe": build_exe_options},
        executables = [Executable("main.py", base=base)])


'''