#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset
import argparse,os,glob
from matplotlib import pyplot as plt
import random
from IPython import embed 
random.seed(123) # make reproducible random number

def readNc(fname):
    fh = Dataset(fname,mode='r')
    n_tim = fh.dimensions['time'].size
    n_lev = fh.dimensions['lev'].size
    CF = np.asarray(fh.variables['CLOUD'][:, :])
    QI = np.asarray(fh.variables['QI'][:,:])
    QL = np.asarray(fh.variables['QL'][:,:])
    QR = np.asarray(fh.variables['QR'][:,:])
    QS = np.asarray(fh.variables['QS'][:,:])
    P = np.asarray(fh.variables['PL'][:,:])
    bscat = np.asarray(fh.variables['backscat'][:,:])
    u = np.asarray(fh.variables['U'][:])
    v = np.asarray(fh.variables['V'][:])
    T = np.asarray(fh.variables['T'][:])
    lat = np.asarray(fh.variables['trjLat'][:])
    lon = np.asarray(fh.variables['trjLon'][:])
    t = np.asarray(fh.variables['time'][:])
    H = np.asarray(fh.variables['H'][:,:])
    fh.close()
        
    return n_tim,n_lev,CF,QI,QL,QR,QS,P/100.,bscat,u,v,T,lat,lon,t,H
def writeNc(fname,nobs,uFlat,vFlat,timezFlat,pFlat,hFlat,angleFlat,latFlat,lonFlat,cfFlat,qiFlat,bFlat):

    fh = Dataset(fname, 'w', format='NETCDF4')
    index  = fh.createDimension('nobs', nobs)

    hours_ = fh.createVariable('hours_in_window', np.float32, ('nobs',))
    hours_[:] = (timezFlat[:]-3.0*3600.)/3600.

    P_ = fh.createVariable('P', np.float32, ('nobs',))
    P_[:] = pFlat[:]

    h_ = fh.createVariable('H', np.float32, ('nobs',))
    h_[:] = hFlat[:]
    
    ang_ = fh.createVariable('ScanAngle', np.float32, ('nobs',))
    ang_[:] = angleFlat[:]

    lat_ = fh.createVariable('trjLat', np.float32, ('nobs',))
    lat_[:] = latFlat[:]

    lon_ = fh.createVariable('trjLon', np.float32, ('nobs',))
    lon_[:] = lonFlat[:]

    u_ = fh.createVariable('U', np.float32, ('nobs',))
    u_[:] = uFlat[:]


    v_ = fh.createVariable('V', np.float32, ('nobs',))
    v_[:] = vFlat[:]
    
    cld_ = fh.createVariable('CLOUD', np.float32, ('nobs',))
    cld_[:] = cfFlat[:]

    qi_ = fh.createVariable('QI', np.float32, ('nobs',))
    qi_[:] = qiFlat[:]

    b_ = fh.createVariable('backscat', np.float32, ('nobs',))
    b_[:] = bFlat[:]

    fh.close()

def filterLidarReturns(CF, QI, QL, QR, QS, bscat, bmin, bmax):
    #colIdx = np.zeros(CF.shape,dtype=bool)
    colIdx = np.zeros(CF.shape,dtype=bool)
    #print(bscat.shape,CF.shape,QI.shape,QL.shape)
 
    for k,cloud_fraction in enumerate(CF):
        if( (cloud_fraction <= 0.05) and (bscat[k]>=bmin) and (bscat[k]<=bmax) and (QI[k]==0) and (QR[k]==0) and (QS[k]==0) ):
            colIdx[k] = True
        elif( (bscat[k]>bmax) and (QR[k]==0) and (QS[k]==0) ):
            colIdx[k] = True
            print('aersol deathblow')
            break
        elif( cloud_fraction>0.0 and QI[k]*1000 <= 0.05 and QI[k]>0.0 and (QL[k]==0) and (QR[k]==0) and (QS[k]==0)):
            colIdx[k] = True
        elif( random.random() < cloud_fraction and cloud_fraction > 0.05 and (QL[k]==0) and (QR[k]==0) and (QS[k]==0) ):
            colIdx[k] = False
            break
    return colIdx          
 

if __name__ == '__main__' :
    parser = argparse.ArgumentParser( description = 'read gz text file,select +/-30 deg scan angle, and output to another file using --in directory and --out directory')
    parser.add_argument('--in', help = 'path to ncdiag', required = True, dest = 'path')
    parser.add_argument('--out', help = 'path to ncdiag', required = True, dest = 'outpath')
    arg = parser.parse_args()
    print("Reading {}".format(arg.path))
    nt,nl,CF,QI,QL,QR,QS,P,bscat,u,v,T,lat,lon,t,H = readNc(arg.path)
    print("Reading done.")
    idx_ret = np.zeros((nt,nl),dtype=bool)
    bmin = 1e-8
    bmax = 100
    print("Filtering.")
    for i in range(nt):
        if(i%1000==0):
            print('index {} of {}'.format(i,nt))
        idx_ret[i,:] = filterLidarReturns(CF[i,:],QI[i,:],QL[i,:],QR[i,:],QS[i,:],bscat[i,:],bmin,bmax)
        if(i%1000==0):
            print('whir')
            print(idx_ret[i,:])
            #print(bscat[i,:])
            print(CF[i,:])
 
    idxOut = np.where(idx_ret)
    lat2 = np.zeros([nt,nl])
    lon2 = np.zeros([nt,nl])
    timez2 = np.zeros([nt,nl])
    angle2 = np.zeros([nt,nl])
    #make flat array with 1d variables extended to match 2d variables with levels
    ang = -30.0
    print("Flattening Data")
    for jjj in range(nt):
        for kkkk in range(nl):
            lat2[jjj,kkkk] = lat[jjj]
            lon2[jjj,kkkk] = lon[jjj]
            timez2[jjj,kkkk] = t[jjj]
            angle2[jjj,kkkk] = ang
        #toggle angle for each provile
        if(ang<0):
            ang = 30.
        else:
            ang = -30.
    latFlat = lat2[idxOut].flatten()
    lonFlat = lon2[idxOut].flatten()
    timezFlat = timez2[idxOut].flatten()
    angleFlat = angle2[idxOut].flatten()
    uFlat = u[idxOut].flatten()
    vFlat = v[idxOut].flatten()
    pFlat = P[idxOut].flatten()
    hFlat = H[idxOut].flatten()
    cfFlat = CF[idxOut].flatten()
    qiFlat = QI[idxOut].flatten()
    bFlat = bscat[idxOut].flatten()
    nobs = latFlat.shape[0]
    print("writting {}".format(arg.outpath))
    writeNc(arg.outpath,nobs,uFlat,vFlat,timezFlat,pFlat,hFlat,angleFlat,latFlat,lonFlat,cfFlat,qiFlat,bFlat)
    print("done!")
