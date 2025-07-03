from setuptools import setup, find_packages

setup(
    name='model_circular_proteins',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'MDAnalysis',
    ],
    entry_points={
        'console_scripts': [
            'build_circle=model_circles.cli:main',
        ],
    },
    author='s-schaef',
    license='BSD 3-Clause',
    description='Tool for building circular protein assemblies.',
)