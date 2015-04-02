from setuptools import setup

__version__ = "0.0.0.8"

setup(name='gaussianemulation',
        version=__version__,
        description='Implementation of gaussian emulation machines for sensitivity analysis..',
        entry_points={},
        url='https://github.com/greenape/gem-module',
        author='Jonathan Gray',
        author_email='j.gray@soton.ac.uk',
        license='MPL',
        packages=['gaussianemulation'],
        include_package_data=True,
        zip_safe=False,
        install_requires=[
        "sympy",
        "mpmath",
    ]
)

with open("gaussianemulation/_version.py", "w") as fp:
        fp.write("__version__ = '%s'\n" % (__version__,))