from setuptools import setup, find_packages

with open('README.md') as fh:
    long_description = fh.read()

setup(
    name='bino',
    version='0.0.1',
    author='Kirill Grishin',
    author_email='grishin@voxastro.org',
    description='Zeropoint calculation for MMIRS and Binospec data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='http://sai.msu.ru',
    license='MIT',
    packages=find_packages(),
    test_suite='test',
    install_requires=['numpy>=1.17.4', 'scipy>=1.3.2', 'matplotlib>=3.0.3', 
                      'pyvo>=1.0', 'argparse>=1.1', 'lmfit>=0.9.15',
                      'astropy>=3.2.3'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: MIT License',
        'Topic :: Education',
        'Programming Language :: Python :: 3',
    ],
    keywords='astrophysics instrumentation',
)