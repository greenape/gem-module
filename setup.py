from numpy.distutils.core import setup, Extension

__version__ = "0.0.1.0"
with open("gaussianemulation/_version.py", "w") as fp:
        fp.write("__version__ = '%s'\n" % (__version__,))

setup(name='gaussianemulation',
        version=__version__,
        description='Implementation of gaussian emulation machines for sensitivity analysis..',
        url='https://github.com/greenape/gem-module',
        author='Jonathan Gray',
        author_email='j.gray@soton.ac.uk',
        license='MPL',
        packages=['gaussianemulation', 'gaussianemulation.gemsa'],
        requires=[
        "sympy",
        "mpmath",
    ],
    ext_modules = [Extension('gaussianemulation.gemsa.gememu', ['gaussianemulation/gemsa/gememu.pyf', 'gaussianemulation/gemsa/gememu.f90'])]
)