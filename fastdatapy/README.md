FAST data reduction - Python

==========


# Intall the package

## download


## intall use pip

	in the folder libfast, open the terminal. 
	
	run: 
	
		python setup.py sdist
	
		cd dist
		
		pip install fastoperations.py-1.0.tar.gz
		
		
	then you can import the libary in jupyter notebook
	
		from fastdatapy.fastoperations import *
		
		
	make sure that all the imported libaries have been intalled, by checking the list '#Import the used libraries' in fastoperations.py
	
	
	


# Pipline in FAST data reduction


## Neptune

extract fits to npy use fast.fitXXYYXYYX2npy(input_path,save_path,finalchan,beam,beam_number)

calibrate, from counts to antenna temperature, use 







This tutorial covers step by step, how to perform a Fast Fourier Transform with Python.

![FFT](https://raw.githubusercontent.com/balzer82/FFT-Python/master/FFT.png)

Including

* How to scale the x- and y-axis in the amplitude spectrum
* Leakage Effect
* Windowing

### [Take a look at the IPython Notebook](http://nbviewer.ipython.org/github/balzer82/FFT-Python/blob/master/FFT-Tutorial.ipynb)

#### Real World Data Example

From

![Vertical Netload Germany 2013](https://raw.githubusercontent.com/balzer82/FFT-Python/master/VerticalGridLoadGermany2013.png)

To

![Periods in NetLoad](https://raw.githubusercontent.com/balzer82/FFT-Python/master/VerticalGridLoadGermany2013-FFT.png)
