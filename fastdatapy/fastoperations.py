#Import the used libraries

import re
import urllib
import sys
import numpy as np
import pandas as pd
from astropy.io import fits
import os
import datetime
import time
from array import array
import matplotlib.pyplot as plt
from pylab import *
import seaborn as sns
from matplotlib.pyplot import MultipleLocator
import scipy.interpolate as spi
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit
import statsmodels.api as sm
import shutil
from pybaselines import Baseline
from pybaselines import utils
from pybaselines.utils import gaussian
from scipy.signal import savgol_filter as sf
from decimal import *
# from matplotlib import cm
# from numba import njit, prange
from scipy.sparse.linalg import spsolve
from scipy import sparse
# from matplotlib.patches import Rectangle

lowess = sm.nonparametric.lowess


#The dates of observation
DATES = np.array(['20210425', '20210426', '20210427', '20210430', '20210501', '20210502',
 '20210503', '20210504', '20210505', '20210506', '20210507', '20210508',
 '20210526'])
 


class fast:




    def fitsXXYYXYYX2npy(self,input_path,save_path,finalchan,beam,beam_number):
        '''
        example:
        input_path = './data/20210526/fits/' #char, the folder path which containts the fits files
        save_path = './npy/20210526/' #char, the folder path which saves the .npy files
        finalchan = 50000 #int, rebin to how many frequency channels you want finally
        beam = '01' # char, which beam: '01' - '19'
        beam_number = 4 #int, total numbers of beam, for computing the mount of fits files in each beam

        '''


        nfile = int(len(os.listdir(input_path))/beam_number)#multi_beam

        finalchan = finalchan

    #     print(math.ceil(nfile/10))


        for n in range (math.ceil(nfile/10)):



            m = 10*n+1

    #         if m == 11:
    #             break

            nlinetotal = 0#update timespec2D.total

            #for i in range (m,m+10):
            for i in range(m,m+10):

                #load fits files, notice the fits' name 
                if i == nfile + 1:
                    break

    #             if i == 11:
    #                 break

                else:
                    hdulist = fits.open(input_path+'Neptune_solarsystracking-M'+beam+'_F_{:>04d}'.format(i)+'.fits')

                print('\rProcess:'+str(i)+'/'+str(nfile),end=' ',flush=True)

                hdu1 = hdulist[1]
                data1 = hdu1.data
                header1 = hdu1.header
                nline = header1['NAXIS2']

    #                 freq = data1['FREQ'][0]
    #                 chwidth = data1['CHAN_BW'][0]
    #                 nchan = data1['NCHAN'][0]
    #                 dt = data1['EXPOSURE'][0]
                #print(freq,chwidth,nchan)

                timespec3D = np.array(data1['DATA'])
                #print(timespec3D.shape)\

    #             #Stokes XX 
                timespec2D_pri = timespec3D[:,:,0]
                #Stokes YY
    #             timespec2D_pri = timespec3D[:,:,1]

    #             return timespec2D_pri


                #remove gain, keep 1050-1450 MHz, need to remove two # following

    #             nchan_pri = int(timespec2D_pri.shape[1] / 10)
                timespec2D = timespec2D_pri#[:,nchan_pri:nchan_pri*9]


                nlinetotal = nlinetotal + nline
                if nlinetotal == nline:
                    timespec2D_total = timespec2D
                else:
                    timespec2D_total = np.concatenate((timespec2D_total,timespec2D),axis = 0)


    #             print(timespec2D_total.shape)


                hdulist.close()#close fits files
                #if i == 1:
                    #break
        #==================================================================================#read and change to array   

      #==================================================================================================
            #combine time channels

            timespec2D = np.transpose(timespec2D_total)
            ntime = 10

            for i in range(0,len(timespec2D[0,:])//10):
                timespec2D_sample = np.mean(timespec2D[:,ntime*i:ntime*(i+1)],axis = 1)
                if i == 0:
                    timespec2D_mean = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_mean = np.concatenate((timespec2D_mean,timespec2D_sample[:,np.newaxis]),axis = 1)

    #         print('time:',timespec2D_mean.shape)



            #interpolate freqency channels for int

            freq_l = len(timespec2D_mean[:,0])

            chan_inter = math.ceil(freq_l/finalchan)*finalchan

            y = np.linspace(0,freq_l,freq_l)

            ynew = np.linspace(0,freq_l,chan_inter)

    #         print(freq_l,chan_inter)

            for t in range(len(timespec2D_mean[0,:])):


                f = spi.interp1d(y,timespec2D_mean[:,t],kind='linear')#(n,m)

                # use linspace so your new range also goes from 0 to 3, with 8 intervals

                timespec2D_mean_interp = f(ynew)

                if t == 0:
                    timespec2D_interp = timespec2D_mean_interp[:,np.newaxis]
                else:
                    timespec2D_interp = np.concatenate((timespec2D_interp,timespec2D_mean_interp[:,np.newaxis]),axis = 1)       


    #         print('interp:',timespec2D_interp.shape)



    #         ##combine frequency chanels================

            nchan = int(timespec2D_interp.shape[0] / finalchan) #1048576channel -> finalchan
    #         print(nchan,timespec2D_total.shape[1] / finalchan)
            for j in range(0,finalchan):
                #print(j)

                timespec2D_sample = np.mean(timespec2D_interp[nchan*j:nchan*(j+1),:],axis = 0)
                if j == 0:
                    timespec2D_com = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_com = np.concatenate((timespec2D_com,timespec2D_sample[:,np.newaxis]),axis = 1)    


            timespec2D_com = timespec2D_com.T


    #         print('frq:',timespec2D_com.shape)

            pth = save_path+beam+'/XX/'
            if not os.path.exists(pth):
                os.makedirs(pth)

            np.save(pth + 'timespec2D_XX_M'+beam+'_'+str(n+1)+'.npy',timespec2D_com)#save as .npy the fast reading file, every ten files#save as .npy the fast reading file, every ten files




        #==================================later for repeating for YY, XY and YX

    #     print(math.ceil(nfile/10))

        for n in range (math.ceil(nfile/10)):



            m = 10*n+1

    #         if m == 11:
    #             break

            nlinetotal = 0#update timespec2D.total

            #for i in range (m,m+10):
            for i in range(m,m+10):

                if i == nfile + 1:
                    break

    #             if i == 11:
    #                 break

                else:
                    hdulist = fits.open(input_path+'Neptune_solarsystracking-M'+beam+'_F_{:>04d}'.format(i)+'.fits')

                print('\rProcess:'+str(i)+'/'+str(nfile),end=' ',flush=True)

                hdu1 = hdulist[1]
                data1 = hdu1.data
                header1 = hdu1.header
                nline = header1['NAXIS2']

    #                 freq = data1['FREQ'][0]
    #                 chwidth = data1['CHAN_BW'][0]
    #                 nchan = data1['NCHAN'][0]
    #                 dt = data1['EXPOSURE'][0]
                #print(freq,chwidth,nchan)

                timespec3D = np.array(data1['DATA'])
                #print(timespec3D.shape)\

    #             #Stokes YY 
    #             timespec2D_pri = timespec3D[:,:,0]
                #Stokes YY
                timespec2D_pri = timespec3D[:,:,1]

    #             return timespec2D_pri


                ##remove gain

    #             nchan_pri = int(timespec2D_pri.shape[1] / 10)
                timespec2D = timespec2D_pri#[:,nchan_pri:nchan_pri*9]
                #print(timespec2D.shape)

                nlinetotal = nlinetotal + nline
                if nlinetotal == nline:
                    timespec2D_total = timespec2D
                else:
                    timespec2D_total = np.concatenate((timespec2D_total,timespec2D),axis = 0)


    #             print(timespec2D_total.shape)


                hdulist.close()#close fits files
                #if i == 1:
                    #break
        #==================================================================================#read and change to array   

      #==================================================================================================
            #combine time channels

            timespec2D = np.transpose(timespec2D_total)
            ntime = 10

            for i in range(0,len(timespec2D[0,:])//10):
                timespec2D_sample = np.mean(timespec2D[:,ntime*i:ntime*(i+1)],axis = 1)
                if i == 0:
                    timespec2D_mean = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_mean = np.concatenate((timespec2D_mean,timespec2D_sample[:,np.newaxis]),axis = 1)

    #         print('time:',timespec2D_mean.shape)



            #interpolate for int

            freq_l = len(timespec2D_mean[:,0])

            chan_inter = math.ceil(freq_l/finalchan)*finalchan

            y = np.linspace(0,freq_l,freq_l)

            ynew = np.linspace(0,freq_l,chan_inter)

    #         print(freq_l,chan_inter)

            for t in range(len(timespec2D_mean[0,:])):


                f = spi.interp1d(y,timespec2D_mean[:,t],kind='linear')#(n,m)

                # use linspace so your new range also goes from 0 to 3, with 8 intervals

                timespec2D_mean_interp = f(ynew)

                if t == 0:
                    timespec2D_interp = timespec2D_mean_interp[:,np.newaxis]
                else:
                    timespec2D_interp = np.concatenate((timespec2D_interp,timespec2D_mean_interp[:,np.newaxis]),axis = 1)       


    #         print('interp:',timespec2D_interp.shape)



    #         ##combine frequency chanels================

            nchan = int(timespec2D_interp.shape[0] / finalchan) #1048576channel -> finalchan
    #         print(nchan,timespec2D_total.shape[1] / finalchan)
            for j in range(0,finalchan):
                #print(j)

                timespec2D_sample = np.mean(timespec2D_interp[nchan*j:nchan*(j+1),:],axis = 0)
                if j == 0:
                    timespec2D_com = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_com = np.concatenate((timespec2D_com,timespec2D_sample[:,np.newaxis]),axis = 1)    

            #np.save('timespec2D_'+str(n+1)+'.npy',timespec2D_sum)# 记得改

            timespec2D_com = timespec2D_com.T


    #         print('frq:',timespec2D_com.shape)

    #         pth = './npy/YY/'
            pth = save_path+beam+'/YY/'

            if not os.path.exists(pth):
                os.makedirs(pth)

            np.save(pth + 'timespec2D_YY_M'+beam+'_'+str(n+1)+'.npy',timespec2D_com)#save as .npy the fast reading file, every ten files#save as .npy the fast reading file, every ten files



        for n in range (math.ceil(nfile/10)):



            m = 10*n+1

    #         if m == 11:
    #             break

            nlinetotal = 0#update timespec2D.total

            #for i in range (m,m+10):
            for i in range(m,m+10):

                if i == nfile + 1:
                    break

    #             if i == 11:
    #                 break

                else:
                    hdulist = fits.open(input_path+'Neptune_solarsystracking-M'+beam+'_F_{:>04d}'.format(i)+'.fits')

                print('\rProcess:'+str(i)+'/'+str(nfile),end=' ',flush=True)

                hdu1 = hdulist[1]
                data1 = hdu1.data
                header1 = hdu1.header
                nline = header1['NAXIS2']

    #                 freq = data1['FREQ'][0]
    #                 chwidth = data1['CHAN_BW'][0]
    #                 nchan = data1['NCHAN'][0]
    #                 dt = data1['EXPOSURE'][0]
                #print(freq,chwidth,nchan)

                timespec3D = np.array(data1['DATA'])
                #print(timespec3D.shape)\

    #             #Stokes XY 
                timespec2D_pri = timespec3D[:,:,2]
                ##remove gain

    #             nchan_pri = int(timespec2D_pri.shape[1] / 10)
                timespec2D = timespec2D_pri#[:,nchan_pri:nchan_pri*9]
                #print(timespec2D.shape)

                nlinetotal = nlinetotal + nline
                if nlinetotal == nline:
                    timespec2D_total = timespec2D
                else:
                    timespec2D_total = np.concatenate((timespec2D_total,timespec2D),axis = 0)


    #             print(timespec2D_total.shape)


                hdulist.close()#close fits files
                #if i == 1:
                    #break
        #==================================================================================#read and change to array   

      #==================================================================================================
            #combine time channels

            timespec2D = np.transpose(timespec2D_total)
            ntime = 10


            for i in range(0,len(timespec2D[0,:])//10):
                timespec2D_sample = np.mean(timespec2D[:,ntime*i:ntime*(i+1)],axis = 1)
                if i == 0:
                    timespec2D_mean = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_mean = np.concatenate((timespec2D_mean,timespec2D_sample[:,np.newaxis]),axis = 1)

    #         print('time:',timespec2D_mean.shape)



            #interpolate for int

            freq_l = len(timespec2D_mean[:,0])

            chan_inter = math.ceil(freq_l/finalchan)*finalchan

            y = np.linspace(0,freq_l,freq_l)

            ynew = np.linspace(0,freq_l,chan_inter)

    #         print(freq_l,chan_inter)

            for t in range(len(timespec2D_mean[0,:])):


                f = spi.interp1d(y,timespec2D_mean[:,t],kind='linear')#(n,m)

                # use linspace so your new range also goes from 0 to 3, with 8 intervals

                timespec2D_mean_interp = f(ynew)

                if t == 0:
                    timespec2D_interp = timespec2D_mean_interp[:,np.newaxis]
                else:
                    timespec2D_interp = np.concatenate((timespec2D_interp,timespec2D_mean_interp[:,np.newaxis]),axis = 1)       


    #         print('interp:',timespec2D_interp.shape)



    #         ##combine frequency chanels================

            nchan = int(timespec2D_interp.shape[0] / finalchan) #1048576channel -> finalchan
    #         print(nchan,timespec2D_total.shape[1] / finalchan)
            for j in range(0,finalchan):
                #print(j)

                timespec2D_sample = np.mean(timespec2D_interp[nchan*j:nchan*(j+1),:],axis = 0)
                if j == 0:
                    timespec2D_com = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_com = np.concatenate((timespec2D_com,timespec2D_sample[:,np.newaxis]),axis = 1)    

            #np.save('timespec2D_'+str(n+1)+'.npy',timespec2D_sum)# 记得改

            timespec2D_com = timespec2D_com.T


    #         print('frq:',timespec2D_com.shape)

            pth = save_path+beam+'/XY/'
            if not os.path.exists(pth):
                os.makedirs(pth)

            np.save(pth + 'timespec2D_XY_M'+beam+'_'+str(n+1)+'.npy',timespec2D_com)#save as .npy the fast reading file, every ten files#save as .npy the fast reading file, every ten files




        #YY

    #     print(math.ceil(nfile/10))

        for n in range (math.ceil(nfile/10)):



            m = 10*n+1

    #         if m == 11:
    #             break

            nlinetotal = 0#update timespec2D.total

            #for i in range (m,m+10):
            for i in range(m,m+10):

                if i == nfile + 1:
                    break

    #             if i == 11:
    #                 break

                else:
                    hdulist = fits.open(input_path+'Neptune_solarsystracking-M'+beam+'_F_{:>04d}'.format(i)+'.fits')

                print('\rProcess:'+str(i)+'/'+str(nfile),end=' ',flush=True)

                hdu1 = hdulist[1]
                data1 = hdu1.data
                header1 = hdu1.header
                nline = header1['NAXIS2']

    #                 freq = data1['FREQ'][0]
    #                 chwidth = data1['CHAN_BW'][0]
    #                 nchan = data1['NCHAN'][0]
    #                 dt = data1['EXPOSURE'][0]
                #print(freq,chwidth,nchan)

                timespec3D = np.array(data1['DATA'])
                #print(timespec3D.shape)\

    #             #Stokes YY 
    #             timespec2D_pri = timespec3D[:,:,0]
                #Stokes YX
                timespec2D_pri = timespec3D[:,:,3]


                ##remove gain

    #             nchan_pri = int(timespec2D_pri.shape[1] / 10)
                timespec2D = timespec2D_pri#[:,nchan_pri:nchan_pri*9]
                #print(timespec2D.shape)

                nlinetotal = nlinetotal + nline
                if nlinetotal == nline:
                    timespec2D_total = timespec2D
                else:
                    timespec2D_total = np.concatenate((timespec2D_total,timespec2D),axis = 0)


    #             print(timespec2D_total.shape)


                hdulist.close()#close fits files
                #if i == 1:
                    #break
        #==================================================================================#read and change to array   

      #==================================================================================================
            #combine time channels

            timespec2D = np.transpose(timespec2D_total)
            ntime = 10


            #print(ntime)
        #     if n == 13:
        #         for i in range(0,125):
        #             timespec2D_sample = np.sum(timespec2D[:,ntime*i:ntime*(i+1)],axis = 1)
        #             if i == 0:
        #                 timespec2D_sum = timespec2D_sample
        #             else:
        #                 timespec2D_sum = np.column_stack((timespec2D_sum,timespec2D_sample))

            for i in range(0,len(timespec2D[0,:])//10):
                timespec2D_sample = np.mean(timespec2D[:,ntime*i:ntime*(i+1)],axis = 1)
                if i == 0:
                    timespec2D_mean = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_mean = np.concatenate((timespec2D_mean,timespec2D_sample[:,np.newaxis]),axis = 1)

    #         print('time:',timespec2D_mean.shape)



            #interpolate for int

            freq_l = len(timespec2D_mean[:,0])

            chan_inter = math.ceil(freq_l/finalchan)*finalchan

            y = np.linspace(0,freq_l,freq_l)

            ynew = np.linspace(0,freq_l,chan_inter)

    #         print(freq_l,chan_inter)

            for t in range(len(timespec2D_mean[0,:])):


                f = spi.interp1d(y,timespec2D_mean[:,t],kind='linear')#(n,m)

                # use linspace so your new range also goes from 0 to 3, with 8 intervals

                timespec2D_mean_interp = f(ynew)

                if t == 0:
                    timespec2D_interp = timespec2D_mean_interp[:,np.newaxis]
                else:
                    timespec2D_interp = np.concatenate((timespec2D_interp,timespec2D_mean_interp[:,np.newaxis]),axis = 1)       


    #         print('interp:',timespec2D_interp.shape)



    #         ##combine frequency chanels================

            nchan = int(timespec2D_interp.shape[0] / finalchan) #1048576channel -> finalchan
    #         print(nchan,timespec2D_total.shape[1] / finalchan)
            for j in range(0,finalchan):
                #print(j)

                timespec2D_sample = np.mean(timespec2D_interp[nchan*j:nchan*(j+1),:],axis = 0)
                if j == 0:
                    timespec2D_com = timespec2D_sample[:,np.newaxis]
                else:
                    timespec2D_com = np.concatenate((timespec2D_com,timespec2D_sample[:,np.newaxis]),axis = 1)    

            #np.save('timespec2D_'+str(n+1)+'.npy',timespec2D_sum)# 记得改

            timespec2D_com = timespec2D_com.T


    #         print('frq:',timespec2D_com.shape)

    #         pth = './npy/YY/'
            pth = save_path+beam+'/YX/'

            if not os.path.exists(pth):
                os.makedirs(pth)

            np.save(pth + 'timespec2D_YX_M'+beam+'_'+str(n+1)+'.npy',timespec2D_com)#save as .npy the fast reading file, every ten files#save as .npy the fast reading file, every ten files



    ### Calibrate   Important!!! 

    #Find Ta

    #new version 2024 Mar.


    ### Calibrate   Important!!! 

    #Find Ta

    #new version 2024 Mar.


    ### Calibrate   Important!!! 

    #Find Ta

    #new version 2024 Mar.


    ### Calibrate   Important!!! 

    #Find Ta

    #new version 2024 Mar.


    def calibrate_count2Ta(self,timespec2D,bandwidth,start_frq,end_frq,save_pth,polarization):
        
        '''
        example:
        timespec2D = np.load('npy/20210526/01/RRL/YY/timespec2D_YY_M01_1.npy') #np.array, 2D array with row as frequency and column as time
        bandwidth = 500 #original bandwidth (1000-1500MHz) of in put array,
        start_frq = 1410 #int, in unit MHz, start freqency from which you want to calibrate 
        end_frq = 1430 #int, in unit MHz, end freqency to which you want to calibrate 
        save_path = 'npy/20210526/01/Ta/' #char, the folder path which saves the atenna temperature .npy files
        polarization = 'YY' # char, which polarization: 'XX','YY','XY' or 'YX'
        
        output:
        figures of noice timespec, noice powerspectrums
        atenna temperature and atenna temperature error .npy
        '''

        Tsys = 25.78 #system temperature

        frq_chan = timespec2D.shape[0]

        #H BW
        Bandwidth = bandwidth #Original bandwidth (1000-1500MHz) of in put array, may change
        start = int((start_frq-1050)/Bandwidth*frq_chan)#start from 1000MHz. If start from 1050, minus 1050
        end = int((end_frq-1050)/Bandwidth*frq_chan)


        timespec2D_sum_frq_t_total = timespec2D

        print(timespec2D_sum_frq_t_total.shape)



    #     count = 0

        for n in range(len(timespec2D_sum_frq_t_total[0,:])//16):
            
            #inject 1-s noice per 16 s, need to change if different

            m = 16*n + 1

            #for i in range(m,m+15):
            ref = timespec2D_sum_frq_t_total[start:end,m-1]-timespec2D_sum_frq_t_total[start:end,m]#discrete


            if m == 1:
                ref_sum = ref[:,np.newaxis]
            else:
                ref_sum = np.concatenate((ref_sum,ref[:,np.newaxis]),axis = 1)

        ref_mean = np.mean(ref_sum,axis = 1)


        plt.figure(figsize=(20,10),dpi = 200,facecolor='white')
        #cmap = sns.cubehelix_palette(reverse = True,as_cmap=True)
        cmap = 'rainbow'
        #sns.heatmap(data = np.log(timespec2D_sum),cmap=sns.cubehelix_palette(as_cmap=True))
        sns.heatmap(data = (ref_sum),cmap = cmap,cbar_kws={'label': 'count'})
        #x_major_locator=MultipleLocator(30)
        y_major_locator=MultipleLocator((end-start)/10)
        ax=plt.gca()
        #ax.xaxis.set_major_locator(x_major_locator)
        ax.yaxis.set_major_locator(y_major_locator)
        #ax.set_xticklabels(['0','0','300','600','900','1200','1500','1800','2100','2400','2700'])
        ax.set_yticklabels([i for i in range(start_frq-int((end_frq-start_frq)/10),end_frq+int((end_frq-start_frq)/10),int((end_frq-start_frq)/10))])

        #plt.xlabel('time(*10s)')
        #plt.ylabel('frequency(*2MHz+1050MHz)')
        plt.xlabel('Time(s)')
        plt.ylabel('Frequency(MHz)')
        plt.title('Noise_Injection:'+polarization)
        #plt.title('Neptune_Temperture'+str(n+1)+':'+str(128*n)+'s-'+str(128*n+49.8)+'s')
        #plt.title('Neptune_cut'+str(n+1)+':'+str(128*n)+'s-'+str(128*n+128)+'s') # 记得改
    #         png_pth = './png/rm_sd/'
    #         if not os.path.exists(pth):
    #             os.makedirs(pth)
    #         plt.savefig(png_pth+'Antena_Temperature:'+polarization[11:13]+'.png') # 记得改
    #         plt.close()
        plt.show()

        frq = np.linspace(start_frq,end_frq,ref_sum.shape[0])
        plt.figure(figsize=(10,3),dpi = 200,facecolor='white')
        plt.plot(frq,np.mean(ref_sum,axis = 1),c='blue',label = 'Original Noice Injection')
        plt.legend()

        baseline_fitter = Baseline(frq, check_finite=False)#creat x axis
        ref_mean_rm_sd = baseline_fitter.arpls(ref_mean, lam=100000.0, diff_order=2, max_iter=50, tol=0.001, weights=None)[0]

        plt.figure(figsize=(10,3),dpi = 200,facecolor='white')
        plt.plot(frq,ref_mean_rm_sd,c='blue',label = 'Noice Injection Smooth once')   
        # print(timespec2D_sum_frq_t_total.shape)        
        # print(ref_mean.shape)
        # Ta_mean = 12.5*(sig - ref)/(2*ref)
        plt.legend()


        ref_mean = baseline_fitter.arpls(-ref_mean_rm_sd, lam=100000.0, diff_order=2, max_iter=50, tol=0.001, weights=None)[0]
        ref_mean = -baseline_fitter.arpls(ref_mean, lam=100000.0, diff_order=2, max_iter=50, tol=0.001, weights=None)[0]

        plt.figure(figsize=(10,3),dpi = 200,facecolor='white')
        plt.plot(frq,ref_mean,c='blue',label = 'Noice Injection Smooth twice') 
        plt.legend()

        #Ta

        #cal every 16 ,2 side

        for n in range(len(timespec2D_sum_frq_t_total[0,:])//16):

            m = 16*n + 1

            for i in range(m,m+15):
                ref = timespec2D_sum_frq_t_total[start:end,m-1]-timespec2D_sum_frq_t_total[start:end,m]#discrete

                if m == 1:
                    sig_ref = timespec2D_sum_frq_t_total[start:end,m-1] - ref#sig except cal at m-1
                else:
                    sig_ref = (timespec2D_sum_frq_t_total[start:end,m] + timespec2D_sum_frq_t_total[start:end,m-2])/2

                Ta_cal = 12.5*sig_ref/ref_mean #temperature at cal
                Ta_cal_err = 12.5*0.02*sig_ref/ref_mean# +  12.5*dC#temperature error at cal


                #===================================================
                #sig = timespec2D_sum_frq_t_total[start:end,i]#i chanel count
                #use single chanel count
                #Ta_temp = 12.5*sig/ref_mean #temperature at no cal

                #===================================================
                #use two chanels count average
                sig = timespec2D_sum_frq_t_total[start:end,i]+timespec2D_sum_frq_t_total[start:end,m-1]#i chanel count
                Ta_temp = 12.5*((sig-ref)/2)/ref_mean
                Ta_temp_err = 12.5*0.02*((sig-ref)/2)/ref_mean# +  12.5*dC#temperature error at cal



                #print(i)
                if i == 1:
                    timespec2D_sum_frq_t_total_Ta = np.column_stack((Ta_cal,Ta_temp))#Ta_total 
                elif i != 1 and i == m:
                    timespec2D_sum_frq_t_total_Ta = np.column_stack((timespec2D_sum_frq_t_total_Ta,Ta_cal,Ta_temp))
                else:
                    timespec2D_sum_frq_t_total_Ta = np.column_stack((timespec2D_sum_frq_t_total_Ta,Ta_temp))

                if i == 1:
                    Ta_err = np.column_stack((Ta_cal_err,Ta_temp_err))#Ta_total 
                elif i != 1 and i == m:
                    Ta_err = np.column_stack((Ta_err,Ta_cal_err,Ta_temp_err))
                else:
                    Ta_err = np.column_stack((Ta_err,Ta_temp_err))



        Ta_pth = save_pth
        if not os.path.exists(Ta_pth):
            os.makedirs(Ta_pth)

    #         print(polarization[11:13])

        np.save(Ta_pth+polarization+'_Ta_'+str(start_frq)+'-'+str(end_frq)+'.npy',timespec2D_sum_frq_t_total_Ta)
        np.save(Ta_pth+polarization+'_Ta_err'+str(start_frq)+'-'+str(end_frq)+'.npy',Ta_err)




        plt.figure(figsize=(20,10),dpi = 200,facecolor='white')
        #cmap = sns.cubehelix_palette(reverse = True,as_cmap=True)
        cmap = 'rainbow'
        #sns.heatmap(data = np.log(timespec2D_sum),cmap=sns.cubehelix_palette(as_cmap=True))
        sns.heatmap(data = (timespec2D_sum_frq_t_total_Ta),cmap = cmap)
        #x_major_locator=MultipleLocator(30)
        y_major_locator=MultipleLocator((end-start)/10)
        ax=plt.gca()
        #ax.xaxis.set_major_locator(x_major_locator)
        ax.yaxis.set_major_locator(y_major_locator)
        #ax.set_xticklabels(['0','0','300','600','900','1200','1500','1800','2100','2400','2700'])
        ax.set_yticklabels([i for i in range(start_frq-int((end_frq-start_frq)/10),end_frq+int((end_frq-start_frq)/10),int((end_frq-start_frq)/10))])
        ax.invert_yaxis()
        ax.set_facecolor("black")
        #plt.xlabel('time(*10s)')
        #plt.ylabel('frequency(*2MHz+1050MHz)')
        plt.xlabel('time(s)')
        plt.ylabel('frequency(MHz)')
        plt.title('Antena_Temperature:'+polarization)
        #plt.title('Neptune_Temperture'+str(n+1)+':'+str(128*n)+'s-'+str(128*n+49.8)+'s')
        #plt.title('Neptune_cut'+str(n+1)+':'+str(128*n)+'s-'+str(128*n+128)+'s') # 记得改
    #         png_pth = './png/rm_sd/'
    #         if not os.path.exists(pth):
    #             os.makedirs(pth)
    #         plt.savefig(png_pth+'Antena_Temperature:'+polarization[11:13]+'.png') # 记得改
    #         plt.close()
        plt.show()













    def dup_rows(self, a, indx, num_dups):
        return np.insert(a,[indx+1]*num_dups,a[indx],axis=0)

    def dup_cols(self, a, indx, num_dups):
        return np.insert(a,[indx+1]*num_dups,a[:,[indx]],axis=1)

    def extract_noise(self,ft2D,gap):

        '''
        Input:
            ft2D: 得到信号的动态谱（二维数组，需要横轴为时间轴，纵轴为频率轴）
            gap: 一个噪声注入周期的长度，比如单位时间是1s, 每 15s 注入一次，持续 1s， 那gap=16

            注：这个算法只适用于 噪声持续时间=单位时间 的情况，假如 ft2D 的时间分辨率是0.1s，而噪声注入 1s, 这个函数还需要改

        Output:
            N_noi_arr_sum: 所有噪声信号的动态谱

        '''

        r = ft2D.shape[1]%gap
        if r == 0:
            n = int(ft2D.shape[1]/gap)
        else:
            n = int(ft2D.shape[1]/gap)+1
        ft2D_ex = ft2D.copy()

        for j in range(0,n):
            N_sig_noi_arr = ft2D[:,j*gap]
            if j == 0:
                N_sig_arr = ft2D[:,j*gap+1]
            else:
                N_sig_arr = (ft2D[:,j*gap+1]+ft2D[:,j*gap-1]) / 2

            N_noi_arr = N_sig_noi_arr-N_sig_arr
            ft2D_ex[:,j*gap] = N_sig_arr

            if j ==0:
                N_noi_arr_sum = N_noi_arr
            else:
                N_noi_arr_sum = np.column_stack((N_noi_arr_sum,N_noi_arr))

        return N_noi_arr_sum,ft2D_ex

    def noise_interp(self,nt,p,arr_2D):
        nf = arr_2D.shape[0]
        nt_n = arr_2D.shape[1]
        x = np.linspace(0,1,2)
        x_new = np.linspace(0,1,p+1)
        y = range(nf)
        for i in range(0,nt_n-1):
            f = interp2d(x,y,arr_2D[:,i:i+2],kind='linear')
            arr_2D_sample_new = f(x_new,y)
            if i == 0:
                arr_2D_new = arr_2D_sample_new[:,:p]
            else:
                arr_2D_new = np.column_stack((arr_2D_new,arr_2D_sample_new[:,:p]))
        arr_sample = arr_2D_new[:,-1]
        for i in range(nt - arr_2D_new.shape[1]):
            arr_2D_new = np.column_stack((arr_2D_new,arr_sample))

        return arr_2D_new

    def index_number(self,li,defaultnumber):
        select = Decimal(str(defaultnumber)) - Decimal(str(li[0]))
        index = 0
        for i in range(1,len(li)-1):
            select2 = Decimal(str(defaultnumber)) - Decimal(str(li[i]))
            if abs(select) > abs(select2):
                select = select2
                index = i
        return [li[index],index]

    def plot_heatmap(self,xaxis,x,y,arr2D,title):

        if xaxis == 'y':
            x_major_locator=MultipleLocator(x)
        elif xaxis == 'n':
            None
        y_major_locator=MultipleLocator(y)
        cmap = 'rainbow'

        plt.figure(figsize=(10,5))
        sns.heatmap(arr2D, cmap=cmap)

        ax=plt.gca()
        if xaxis == 'y':
            ax.xaxis.set_major_locator(x_major_locator)
        elif xaxis == 'n':
            None

        ax.yaxis.set_major_locator(y_major_locator)

        if xaxis == 'y':
            ax.set_xticklabels(['0','0','500','1000','1500','2000','2500'])
            plt.xlabel('Time [s]')
        elif xaxis == 'n':
            plt.xlabel('Number of noise injection')

        ax.set_yticklabels(['1410','1410','1414','1418','1422','1426','1430'])
        plt.ylabel('Frequency [MHz]')
        plt.title(title)
        plt.show()

    def bin1D(self,arr,newlen):
        len = arr.shape[0]
        list_new = []
        if len % newlen == 0:
            bin = int(len / newlen)
            for i in range(newlen):
                list_new.append(np.mean(arr[i*bin:(i+1)*bin]))
        else:
            len_in = (int(len/newlen)+1)*newlen
            bin = int(len_in / newlen)
            f = interp1d(range(len),arr)
            arr_in = f(np.linspace(0,len-1,len_in))
            for i in range(newlen):
                list_new.append(np.mean(arr_in[i*bin:(i+1)*bin]))

        return np.asarray(list_new)

    def bin2D(self,ax,arr2D,nbin):
        nr, nc = arr2D.shape
        x = range(nc)
        y = range(nr)
        if ax == 0:
            nchan = int(nr/nbin)+1
            nr_new = (nchan)*nbin
            y_new = np.linspace(0,nr-1,nr_new)

            f = interp2d(x,y,arr2D)
            arr2D_new = f(x,y_new)
            # for i in range(nc):
            #     f = interp1d(x,arr2D[:,i])
            #     arr1D_new = f(x_new)
            #     if i == 0:
            #         arr2D_new = arr1D_new
            #     else:
            #         arr2D_new = np.column_stack((arr2D_new,arr1D_new))

            for j in range(nbin):
                sample = np.mean(arr2D_new[nchan*j:nchan*(j+1),:],axis=0)
                if j == 0:
                    arr2D_bin = sample
                else:
                    arr2D_bin = np.row_stack((arr2D_bin,sample))

        if ax == 1:
            nchan = int(nc/nbin)+1
            nc_new = (nchan)*nbin
            x_new = np.linspace(0,nc-1,nc_new)

            f = interp2d(x,y,arr2D)
            arr2D_new = f(x_new,y)
            # for i in range(nr):
            #     f = interp1d(x,arr2D[i,:])
            #     arr1D_new = f(x_new)
            #     if i == 0:
            #         arr2D_new = arr1D_new
            #     else:
            #         arr2D_new = np.row_stack((arr2D_new,arr1D_new))

            for j in range(nbin):
                sample = np.mean(arr2D_new[:,nchan*j:nchan*(j+1)],axis=1)
                if j == 0:
                    arr2D_bin = sample
                else:
                    arr2D_bin = np.column_stack((arr2D_bin,sample))


        return arr2D_bin

    def delta_corr_2D(self,ci_2D,cr_2D,delta_2D):

        '''
        对arctan得到的结果delta_2D进行修正

        Input:
            cr_2D: CR 的二维数组
            ci_2D: CI 的二维数组
            delta_2D: CI/CR 可得tan\delta，做反三角函数可得 delta，但这个 delta 的象限需要进一步确定

        Output:
            更新好象限位置的 delta

        '''
        # [pi/2,pi]
        mask_sin_p = ci_2D > 0
        # print(mask_sin_p)
        mask_cos_n = cr_2D < 0
        mask_sin_p_cos_n = np.logical_and(mask_sin_p,mask_cos_n)
        delta_2D[mask_sin_p_cos_n] = delta_2D[mask_sin_p_cos_n] + np.pi

        # [-pi,-pi/2]
        mask_sin_n = ci_2D < 0
        mask_sin_n_cos_n = np.logical_and(mask_sin_n,mask_cos_n)
        delta_2D[mask_sin_n_cos_n] = delta_2D[mask_sin_n_cos_n] - np.pi
        return delta_2D 

    def ArPLS(self, y, lam, ratio=0.05, itermax=10):
        '''
        copy from https://irfpy.irf.se/projects/ica/_modules/irfpy/ica/baseline.html

        Baseline correction using asymmetrically
        reweighted penalized least squares smoothing
        Sung-June Baek, Aaron Park, Young-Jin Ahna and Jaebum Choo,
        Analyst, 2015, 140, 250 (2015)

        Inputs:
            y:
                input data (i.e. SED curve)
            lam:
                parameter that can be adjusted by user. The larger lambda is,
                the smoother the resulting background, z
            ratio:
                wheighting deviations: 0 < ratio < 1, smaller values allow less negative values
            itermax:
                number of iterations to perform initial:10
        Output:
            the fitted background vector
        '''

        N = len(y)
        #  D = sparse.csc_matrix(np.diff(np.eye(N), 2))
        D = sparse.eye(N, format='csc')
        D = D[1:] - D[:-1]  # numpy.diff( ,2) does not work with sparse matrix. This is a workaround.
        D = D[1:] - D[:-1]

        D = D.T
        w = np.ones(N)
        for i in range(itermax):
            W = sparse.diags(w, 0, shape=(N, N))
            Z = W + lam * D.dot(D.T)
            z = spsolve(Z, w * y)
            d = y - z
            dn = d[d < 0]
            m = np.mean(dn)
            s = np.std(dn)
            wt = 1. / (1 + np.exp(2 * (d - (2 * s - m)) / s))
            if np.linalg.norm(w - wt) / np.linalg.norm(w) < ratio:
                break
            w = wt

        return z

    def WhittakerSmooth(self, x,w,lambda_,differences=1):
        '''
        Penalized least squares algorithm for background fitting

        input
            x: input data (i.e. chromatogram of spectrum)
            w: binary masks (value of the mask is zero if a point belongs to peaks and one otherwise)
            lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background
            differences: integer indicating the order of the difference of penalties

        output
            the fitted background vector
        '''
        X=np.matrix(x)
        m=X.size
        E= sparse.eye(m,format='csc')
        for i in range(differences):
            E=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
        W=sparse.diags(w,0,shape=(m,m))
        A=sparse.csc_matrix(W+(lambda_*E.T*E))
        B=sparse.csc_matrix(W*X.T)
        background=spsolve(A,B)
        return np.array(background)

    def airPLS(self, x, lambda_=100, porder=1, itermax=15):
        '''
        Adaptive iteratively reweighted penalized least squares for baseline fitting

        input
            x: input data (i.e. chromatogram of spectrum)
            lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background, z
            porder: adaptive iteratively reweighted penalized least squares for baseline fitting

        output
            the fitted background vector
        '''
        m=x.shape[0]
        w=np.ones(m)
        for i in range(1,itermax+1):
            z=WhittakerSmooth(x,w,lambda_, porder)
            d=x-z
            dssn=np.abs(d[d<0].sum())
            if(dssn<0.001*(abs(x)).sum() or i==itermax):
                if(i==itermax): print('WARING max iteration reached!')
                break
            w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
            w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
            w[0]=np.exp(i*(d[d<0]).max()/dssn) 
            w[-1]=w[0]
        return z

    def median_filter_1D(self,arr,wl):
        '''
        wl: window length, must be odd
        arr: 1D array 
        '''
        if wl%2 == 0:
            print('Window length must be odd!')
            return None
        else:
            # add
            nf = arr.shape[0]
            nf_new = nf + (wl-1)
            d1 = int((wl-1)/2)
            d2 = int((wl+1)/2)
            arr_new = np.zeros((nf_new))
            arr_new[:d1] = np.ones((d1))*np.mean(arr[:d2])
            arr_new[d1:(-1)*d1] = arr
            arr_new[(-1)*d1:] = np.ones((d1))*np.mean(arr[(-1)*d2:])
            # print(arr_new)

            median_list = list()
            for i in range(nf):
                median = np.median(arr_new[i:(i+wl)])
                median_list.append(median)
            median_arr = np.asarray(median_list)
            return median_arr

    def smooth_fs_2D(self,arr_exn,method,l):

        '''
        For standing waves or RFI removal.
        You could choose the smoothing method.
        Data must be: axis_0--Frequency, axis_1--Time
        '''

        nt = arr_exn.shape[1]
        arr_exn_sm = arr_exn.copy()

        if method == 'ArPLS':
            for i in range(nt):
                arr_exn_sm[:,i] = ArPLS(arr_exn[:,i],l)
        elif method == 'median':
            for i in range(nt):
                arr_exn_sm[:,i] = median_filter_1D(arr_exn[:,i],l)
        else:
            print('method NAME is wrongggggg')

        return arr_exn_sm

    def remove_baseline_1D(self,raw_arr,on_off_swiches,quasi_cycle_length,mini_gap):

        '''
        output: arg_begin, arg_end, baseline_arr, ebl_arr
        '''

        # extract the points which are on the baseline
        baseline_list = []
        idx_baseline_list = []
        for i in range(on_off_swiches):
            if i == 0:
                idx_cycle_start = i*quasi_cycle_length
            else:
                idx_cycle_start = max([idx_baseline_list[i-1]+mini_gap,i*quasi_cycle_length])
            idx_cycle_end = (i+1)*quasi_cycle_length
            print(idx_cycle_start,idx_cycle_end)
            cycle = raw_arr[idx_cycle_start:idx_cycle_end]
            baseline_point = np.min(cycle)
            idx_baseline_point = np.argmin(cycle)+idx_cycle_start

            baseline_list.append(baseline_point)
            idx_baseline_list.append(idx_baseline_point)

        # the first and last points redefine the length of the baseline-romoved data
        idx_baseline_start = idx_baseline_list[0]
        idx_baseline_end = idx_baseline_list[-1]+1

        # interpolation for the baseline
        f = interp1d(idx_baseline_list,baseline_list)
        idx_baseline_list_interp = np.arange(idx_baseline_start,idx_baseline_end)
        baseline = f(idx_baseline_list_interp)

        # remove the baseline
        new_arr = raw_arr[idx_baseline_start:idx_baseline_end]-baseline

        return idx_baseline_list,baseline,new_arr

    def remove_baseline_1D_known(self,raw_arr,index_baseline_arr_1D):

        baseline_list = []
        for i in index_baseline_arr_1D:
            baseline_list.append(raw_arr[i])

        # the first and last points redefine the length of the baseline-romoved data
        baseline_start_index = index_baseline_arr_1D[0]
        baseline_end_index = index_baseline_arr_1D[-1]+1

        # interpolation for the baseline
        f = interp1d(index_baseline_arr_1D,baseline_list)
        index_baseline_list_interp = np.arange(baseline_start_index,baseline_end_index)
        baseline = f(index_baseline_list_interp)

        # remove the baseline
        new_arr = raw_arr[baseline_start_index:baseline_end_index]-baseline
        return baseline,new_arr

    def remove_baseline_2D(self,arr_exn,on_off_swiches,quasi_cycle_length,mini_gap):

        baseline_start_index_list = []
        baseline_end_index_list = []
        new_arr_1D_list = []
        nf = arr_exn.shape[0]

        for j in range(nf):
            index_baseline_list,baseline,new_arr_1D = remove_baseline_1D(arr_exn[j,:],on_off_swiches,quasi_cycle_length,mini_gap)
            # save the baseline points index in an array for Stokes L and delta_s
            if j == 0:
                baseline_point_index_arr = np.asarray(index_baseline_list)
            else:
                baseline_point_index_arr = np.row_stack((baseline_point_index_arr,np.asarray(index_baseline_list)))

            # the first and last points redefine the length of the baseline-romoved data
            baseline_start_index = index_baseline_list[0]
            baseline_end_index = index_baseline_list[-1]+1
            # and for a 2d-array, the length of every channel should be the same
            baseline_start_index_list.append(baseline_start_index)
            baseline_end_index_list.append(baseline_end_index)
            new_arr_1D_list.append(new_arr_1D)

        # find the shared length
        new_arr_exn_start_index = max(baseline_start_index_list)
        new_arr_exn_end_index = min(baseline_end_index_list)
        # and cutting for the harmony
        for k in range(nf):
            d_start = new_arr_exn_start_index - baseline_start_index_list[k]
            d_end = new_arr_exn_end_index - baseline_end_index_list[k]
            if d_end == 0:
                new_arr_1D = new_arr_1D_list[k][d_start:]
            else:
                new_arr_1D = new_arr_1D_list[k][d_start:d_end]

            if k == 0:
                new_arr_exn = new_arr_1D
            else:
                new_arr_exn = np.row_stack((new_arr_exn,new_arr_1D))

        return baseline_point_index_arr,new_arr_exn_start_index,new_arr_exn_end_index,new_arr_exn

    def remove_baseline_2D_known(self,arr_exn,baseline_point_index_arr):
        baseline_start_index_list = []
        baseline_end_index_list = []
        new_arr_1D_list = []
        nf = arr_exn.shape[0]

        for j in range(nf):
            index_baseline_arr_1D = baseline_point_index_arr[j,:]
            baseline,new_arr_1D = remove_baseline_1D_known(arr_exn[j,:],index_baseline_arr_1D)
            # the first and last points redefine the length of the baseline-romoved data
            baseline_start_index = index_baseline_arr_1D[0]
            baseline_end_index = index_baseline_arr_1D[-1]+1
            # and for a 2d-array, the length of every channel should be the same
            baseline_start_index_list.append(baseline_start_index)
            baseline_end_index_list.append(baseline_end_index)
            new_arr_1D_list.append(new_arr_1D)

        # find the shared length
        new_arr_exn_start_index = max(baseline_start_index_list)
        new_arr_exn_end_index = min(baseline_end_index_list)
        # and cutting for the harmony
        for k in range(nf):
            d_start = new_arr_exn_start_index - baseline_start_index_list[k]
            d_end = new_arr_exn_end_index - baseline_end_index_list[k]
            if d_end == 0:
                new_arr_1D = new_arr_1D_list[k][d_start:]
            else:
                new_arr_1D = new_arr_1D_list[k][d_start:d_end]

            if k == 0:
                new_arr_exn = new_arr_1D
            else:
                new_arr_exn = np.row_stack((new_arr_exn,new_arr_1D))

        return new_arr_exn_start_index,new_arr_exn_end_index,new_arr_exn
