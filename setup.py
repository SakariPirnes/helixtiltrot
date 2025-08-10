import setuptools


with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='helixside',
    version='0.0.2',
    author='Sakari Pirnes',
    author_email='sakari.pirnes@helsinki.fi',
    description='nice stuff',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    packages=['helixside']

)
