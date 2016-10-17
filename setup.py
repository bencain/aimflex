from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='aimflex',
      version='0.1',
      description='Analytic Image Modeling with Flexion',
      url='http://github.com/bencain/aimflex',
      author='Benjamin M Cain',
      author_email='benjamin.m.cain@gmail.com',
      license='MIT',
      packages=['aimflex'],
      install_requires=[
      	'astropy',
      	'numpy',
      	'emcee'
      ],
      zip_safe=False)
