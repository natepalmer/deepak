from setuptools import setup

setup(name='deepak',
      version='0.1',
      description='Analysis suite for deep mutational scan data',
      url='http://github.com/natepalmer/deepak',
      author='Nathan Palmer',
      author_email='ndpalmer@ucsd.edu',
      license='MIT',
      packages=['deepak'],
      install_requires=[
          'numpy', 'scipy', 'matplotlib', 'pandas', 'biopython'
      ],
      entry_points={
            'console_scripts': ['deepak-classify_wrapper=deepak.main:classify',
                                'deepak-quantify=deepak.main:quantify'],
      },
      zip_safe=False)
