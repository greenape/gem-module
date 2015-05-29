from numpy.distutils.core import setup, Extension

__version__ = "0.0.0.9"

setup(name='gaussianemulation',
        version=__version__,
        description='Implementation of gaussian emulation machines for sensitivity analysis..',
        url='https://github.com/greenape/gem-module',
        author='Jonathan Gray',
        author_email='j.gray@soton.ac.uk',
        license='MPL',
        packages=['gaussianemulation'],
        requires=[
        "sympy",
        "mpmath",
    ],
    ext_modules = [Extension('gaussianemulation.gememu', ['gaussianemulation/gemsa/gememu.pyf', 'gaussianemulation/gemsa/gememu.f90'])]
)

with open("gaussianemulation/_version.py", "w") as fp:
        fp.write("__version__ = '%s'\n" % (__version__,))