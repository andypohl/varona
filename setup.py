import pathlib

import setuptools


# Function to read the version from __version__.py
def read_version():
    version = {}
    version_file = pathlib.Path(__file__).parent / "src" / "verona" / "__version__.py"
    with version_file.open() as f:
        exec(f.read(), version)
    return version["__version__"]


setuptools.setup(version=read_version())
