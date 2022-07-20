try:
    from skbuild import setup
except ImportError:
    raise ImportError('scikit-build is required for installing from source')

setup(
    name="spherical-rtree",
    version="0.1.0",
    description="Spherical RTree for point-line relationship queries",
    author='Jacob Zhong',
    author_email='cmpute@gmail.com',
    url='https://github.com/cmpute/spherical-rtree',
    download_url='https://github.com/cmpute/spherical-rtree/archive/master.zip',
    license='MIT',
    packages=['sphertree'],
    install_requires=['numpy>=1.11'],
    setup_requires=['cython>=0.29', 'scikit-build'],
    extras_require={'test': ['pytest']},
    package_data={'spherical_rtree': []},
    classifiers=[
        'Programming Language :: C++',
        'Programming Language :: Cython',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Development Status :: 2 - Pre-Alpha',
        'Topic :: Scientific/Engineering'
    ],
    keywords=['rtree', 'spherical', 'spatial index', 'cython'],
)