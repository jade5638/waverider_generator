from setuptools import setup, find_packages

VERSION = '1.0.0'
DESCRIPTION = 'Hypersonic Waverider Generator'
LONG_DESCRIPTION = 'Package which allows a user to generate and export parametric waverider geometries'

setup(
    name="waverider_generator",
    version=VERSION,
    author="Jade Nassif",
    author_email="<jade.nassif2002@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url='https://github.com/jade5638/waverider_generator',
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['numpy','cadquery','scipy','matplotlib'],
    keywords=['python', 'waverider', 'CAD', 'hypersonic','aerospace-engineering']
)
