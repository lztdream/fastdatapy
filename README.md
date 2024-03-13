# fastdatapy
FAST data reduction - Python


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

calibrate, from counts to antenna temperature
