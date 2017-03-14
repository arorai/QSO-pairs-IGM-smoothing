import numpy as np
import os
'''
The inputs required for building an emulator are: 

1-A parameter grid, in the form of an ASCII table

2-The statistic you need to emulate, in the form 
of a n_m x n_p matrix, where n_p is the number 
of points of the emulated function (e.g. the number
of k-bins of the power spectrum) and n_m is the 
number of models in the grid. The ordering of the 
models must be consistent with the provided 
parameter grid.

3-you need to specify the smoothing lengths for 
the gaussian process: these parameters express
the degree of smoothness of the statistic along
each direction in parameter space

4-you may want to transform the parameter
space if you think that the emulated statistic
behave better as a function of different physical
parameters than those specified in the grid

5-specify a name which will be used to store files

###################################################
OUTPUT
During the elaboration the code will store:

-The parameter grid
-The principal components of the statistic + mean and variance
-The weights of each component for each point in the gread

-The correlation matrix used by the Gaussian process
-The Gaussian process smoothing lengths you chose

This files are what the emulator requires to work properly. 
Then you can create your emulator calling the function

Myemulator = create_emulator('stat_name')

where 'stat_name' is the name you assigned at point 5 above.
Myemulator will be a function which interpolates your statistic
in parameter space 

statistics_in_x = Myemulator( x , n_PC )

where x is the desired point in parameter space, and n_PC
the desired number of principal component

###################################################
PACKAGE REQUIRED
Beside numpy, os and pickle, which are probably already
included in your python version, you need to install infpy which 
contains the gaussian processes routines 



'''
def emulator_builder(par_grid,stat_grid, name,em_path, smooth_lengths):
	'''
	create folder, if does not exist already
	'''
	path=em_path+name+'_emulator/'
	if not os.path.exists(path):
		os.makedirs(path)
	file=open(path+name+'_parameter_grid.txt','w')
	pickle.dump(par_grid,file)
	file.close()
	'''
	first calculate and store PCAs
	'''
	PC,ww,mean_stat,var_stat=PCA_analysis(stat_grid) #function incomplete
	
	file=open(path+name+'_PC.txt','w')
	pickle.dump(PC,file)
	file.close()
        file=open(path+name+'_ww.txt','w')
        pickle.dump(ww,file)
        file.close()
        file=open(path+name+'_mean.txt','w')
        pickle.dump(mean_stat,file)
        file.close()
        file=open(path+name+'_var.txt','w')
        pickle.dump(var_stat,file)
        file.close()
	'''	
	calculate correlation matrix
	'''
	npcmax=min([len(PC),50])
	Rwstar=Rwstar_matrix(par_grid,ww,smooth_lengths,npcmax) 
 
        file=open(path+name+'_Rwstar.txt','w')
        pickle.dump(Rwstar,file)
        file.close()
        file=open(path+name+'_smooth_lengths.txt','w')
        pickle.dump(smooth_lengths,file)
        file.close()

	return


def create_emulator(name,em_path):
	#'em_path' is the folder where the emulators are stored
	path=em_path+name+'_emulator/'

	file=open(path+name+'_parameter_grid.txt','r')
	par_grid=pickle.load(file)
	file.close()

	file=open(path+name+'_PC.txt','r')
	PC=pickle.load(file)
	file.close()

        file=open(path+name+'_ww.txt','r')
        ww=pickle.load(file)
        file.close()

        file=open(path+name+'_mean.txt','r')
        mean_stat=pickle.load(file)
        file.close()

        file=open(path+name+'_var.txt','r')
        var_stat=pickle.load(file)
        file.close()

	file=open(path+name+'_Rwstar.txt','r')
        Rwstar=pickle.load(file)
        file.close()

        file=open(path+name+'_smooth_lengths.txt','r')
        smooth_lengths=pickle.load(file)
        file.close()

	predict_weights=build_gp_predictor(smooth_lengths,Rwstar,par_grid)
	emulator=build_emulator(mean_stat,var_stat,PC,predict_weights)
	return emulator
