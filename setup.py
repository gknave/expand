from setuptools import setup
setup(name='expand_nd',
    version='0.1',
	description='Scratch package for methods to expand an n-dimensional ODE.',
	url='https://github.com/gknave/expand_nd',
	author='Gary Nave',
	author_email='gknave@vt.edu',
	packages=['expand_nd'],
	install_requires=[
	  'numpy',
	  'sympy',
	],
	zip_safe=False)