import numpy as np
import pickle
import infpy
from infpy.gp import GaussianProcess
from infpy.gp import SquaredExponentialKernel as SE
from infpy.gp import NeuralNetworkKernel as NN
from infpy.gp import noise_kernel as noise



def SoftenedSquaredExponentialKernel(smooth_lengths,soft_l):
	def kernel(x1,x2,identical=False):
		dist=(x1-x2)/smooth_lengths
		d=max([np.dot(dist,dist),1.2])
		if identical: return 1
		else:	return exp(-d**2/2)
	return kernel

def PCA_analysis(stat_grid):
	if np.ndim(stat_grid)==1:
		meanst=np.mean(stat_grid)
		varst=var(stat_grid)
		renorm=(stat_grid-meanst)/sqrt(varst)
		return [1],renorm,meanst,varst
	nbins=len(stat_grid[0,:])
	nm=len(stat_grid[:,0])
	npc=min(nm,nbins)
	meangrid=np.array([mean(stat_grid[:,i]) for i in xrange(nbins)])
	vargrid=np.array([var(stat_grid[:,i]) for i in xrange(nbins)])
	#WATCH OUT! IF THE VARIANCE IS 0 THE RENORMALIZATION
	#BELOW WILL MESS UP EVERYTHING!!!
	vargrid[where(vargrid == 0)]=1.0 #this should fix it?
	normgrid=np.array([(stat_grid[:,i]-meangrid[i])/sqrt(vargrid[i]) for i in xrange(nbins)])
	newgrid=reshape(normgrid,(nbins,nm))
	#svd and pca decomposition
	U,D,V=svd(newgrid)
	PC=np.array([(U[:,i]*D[i]) for i in xrange(npc)])#these are the principal functions
	ww=V #these are the relative weights
	return PC,ww,meangrid,vargrid
###################################################
def Rwstar_matrix(X,ww,Smooth_lengths,npc):

	ndim=len(Smooth_lengths)
	npts=len(X[:,0])
	Rw_stars = matrix(zeros((npts,npc)))
	for i in xrange(npc):
		if npc==1:
			Y=ww
		else:
			Y=ww[i,:npts]
		kernel=SE(params=Smooth_lengths,dimensions=ndim)+noise([0.15])	
		#kernel=SoftenedSquaredExponentialKernel(Smooth_lengths,1)
		gpnow=GaussianProcess(X,Y,kernel)
		Rw_stars[:,i]=np.array(gpnow._alpha)
	return Rw_stars
#############################
def build_gp_predictor(Smooth_lengths,Rwstar,X):
	npts=len(X[:,0])
	def predict_weights(x,npc):
	        #wnow=np.array([3-log10(x[0]),1-0.5*x[1],x[2]]  )
                diff= (X-x[np.newaxis,:])/Smooth_lengths[np.newaxis,:]              
                k_new=np.exp(-(diff*diff).sum(-1)/2.)
                new_weights= np.array(np.dot(k_new,Rwstar))[0,:npc]
		
	        #new_weights=zeros(npc,float)
	        #for j in xrange(npc):
		#	diff=np.array([(X[i,:]-x)/Smooth_lengths for i in xrange(npts)])
	        #        newk=[exp(-np.dot(diff[i],diff[i])/2) for i in xrange(npts)]
	        #        k_new=np.array(newk)
		#	new_weights[j]= float(np.dot(k_new,Rwstar[:,j]))
		
	        return new_weights
	return predict_weights
#############################
def build_emulator(mean,var,pca,predict_new_weights):
	def emulator(x,npc):
		em=0
	        new_weights=predict_new_weights(x,npc)
		#em=np.sum(np.sqrt(var)[np.newaxis,:]*pca[0:npc,:]*new_weights[:,np.newaxis],axis=0)
	        for i in xrange(npc):
        	    em+=sqrt(var)*pca[i]*new_weights[i]
        	return em+mean		
	return emulator



