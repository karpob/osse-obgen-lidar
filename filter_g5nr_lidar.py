#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset
import argparse,os,glob
from matplotlib import pyplot as plt
import random
from IPython import embed 
random.seed(123) # make reproducible random number
#units from https://www.nwpsaf.eu/site/download/documentation/rtm/docs_rttov12/rttov_gas_cloud_aerosol_units.pdf
Mair = 28.9648
Mh2o = 18.01528
epsilon = Mair/Mh2o
Mdry = 28.9644
R = 8.3144598 
Rd = R/Mdry
def main(inArg,outArg):
    nt,nl,CF,QI,QL,QR,QS,QV,P,bscat,u,v,T,lat,lon,t,H = readNc(inArg)
    print("Reading done.")
    idx_ret = np.zeros((nt,nl),dtype=float)
    bmin = 1e-6
    bmax = 1e2
    print("Filtering.")
    
    QII = mass2VolumeDensity(QI,P,T, QV)
    QLL = mass2VolumeDensity(QL,P,T,QV)
    for i in range(nt):
        if(i%1000==0):
            print('index {} of {}'.format(i,nt))
        idx_ret[i,:] = filterLidarReturns(CF[:,:],QII[i,:],QLL[i,:],QR[i,:],QS[i,:],QV[i,:], P[i,:], T[i,:],H[i,:],bscat[:,:],bmin,bmax,i)
    print(idx_ret.min(),idx_ret.max())
    idxOut = np.where(idx_ret>0)
    lat2 = np.zeros([nt,nl])
    lon2 = np.zeros([nt,nl])
    timez2 = np.zeros([nt,nl])
    levz2 = np.zeros([nt,nl])
    levz = np.arange(1,73)
    #make flat array with 1d variables extended to match 2d variables with levels
    print("Flattening Data")
    for jjj in range(nt):
        for kkkk in range(nl):
            lat2[jjj,kkkk] = lat[jjj]
            lon2[jjj,kkkk] = lon[jjj]
            timez2[jjj,kkkk] = t[jjj]
            levz2[jjj,kkkk] = levz[kkkk]
    latFlat = lat2[idxOut].flatten()
    lonFlat = lon2[idxOut].flatten()
    timezFlat = timez2[idxOut].flatten()
    uFlat = u[idxOut].flatten()
    vFlat = v[idxOut].flatten()
    pFlat = P[idxOut].flatten()
    hFlat = H[idxOut].flatten()
    cfFlat = CF[idxOut].flatten()
    qiFlat = QI[idxOut].flatten()
    bFlat = bscat[idxOut].flatten()
    levzFlat = levz2[idxOut].flatten()
    idxFlat = idx_ret[idxOut].flatten()
    nobs = latFlat.shape[0]
    print("writting {}".format(outArg))
    writeNc(outArg,nobs,uFlat,vFlat,timezFlat,pFlat,hFlat,latFlat,lonFlat,cfFlat,qiFlat,bFlat,levzFlat,idxFlat)
    print("done!")


def moistR(qv):
    mR = Rd*(1.0+ ( ( (1-epsilon)/epsilon )*qv))
    return mR
    
def mass2VolumeDensity(qil,p,t,qv):
    """
    kgkg of cloud from kg/kg to g/m3
    """
    
    RT = moistR(qv) * t
    WC =  qil*(p/RT)

    return WC


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
    QV = np.asarray(fh.variables['QV'][:,:])
    fh.close()
        
    return n_tim,n_lev,CF,QI,QL,QR,QS,QV,P,bscat,u,v,T,lat,lon,t,H
def writeNc(fname,nobs,uFlat,vFlat,timezFlat,pFlat,hFlat,latFlat,lonFlat,cfFlat,qiFlat,bFlat,levzFlat,idxFlat):

    fh = Dataset(fname, 'w', format='NETCDF4')
    index  = fh.createDimension('nobs', nobs)

    hours_ = fh.createVariable('hours_in_window', np.float32, ('nobs',))
    hours_[:] = (timezFlat[:]-3.0*3600.)/3600.

    P_ = fh.createVariable('P', np.float32, ('nobs',))
    P_[:] = pFlat[:]

    h_ = fh.createVariable('H', np.float32, ('nobs',))
    h_[:] = hFlat[:]
    
    lat_ = fh.createVariable('trjLat', np.float32, ('nobs',))
    lat_[:] = latFlat[:]

    lon_ = fh.createVariable('trjLon', np.float32, ('nobs',))
    lon_[:] = lonFlat[:]

    print('nobs',nobs,levzFlat.shape)
    levz_ = fh.createVariable('Level', np.float32, ('nobs',))
    levz_[:] = levzFlat[:]

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

    iii_ = fh.createVariable('returnType',np.int,('nobs',))
    iii_[:] = idxFlat[:]
    fh.close()

def cloudFracPdf(cloud_fraction):
   
    if(cloud_fraction < 0.2):
        pCloud = cloud_fraction
    elif(cloud_fraction<0.4):
        pCloud = 0.35
    elif(cloud_fraction <0.6):
        pCloud = 0.55
    elif(cloud_fraction<0.8):
        pCloud = 0.75
    elif(cloud_fraction<0.9):
        pCloud = 0.85
    elif(cloud_fraction<0.95):
        pCloud = 0.90
    elif(cloud_fraction>0.95 and cloud_fraction<1.0):
        pCloud = 0.95
    elif( np.isclose(cloud_fraction,1.0) ):
        pCloud = 0.99
    
    pCloud=cloud_fraction
    return pCloud


def filterLidarReturns(CF, QI, QL, QR, QS, QV, P, T,H, bbscat, bmin, bmax, ix):
    
    colIdx = np.zeros(CF.shape[1],dtype=int)
    
    ixMax = CF.shape[0]
    izMax = CF.shape[1]
    for k,cloud_fraction in enumerate(CF[ix,:]):
        pCloud = cloudFracPdf(cloud_fraction)
        qivol = QI[k]
        qlvol = QL[k]
        bscat = bbscat[ix,:]
       

        #if clear and sufficient aerosol set to True
        if( (cloud_fraction <= 0.05) and (bscat[k]>bmin) and (bscat[k]<bmax) and np.isclose(qivol,0) and np.isclose(qlvol,0) and np.isclose(QR[k],0) and np.isclose(QS[k],0) ):
            colIdx[k] = 1
        #if we have a thin cirrus treat it like an aerosol in clear conditions
        elif( cloud_fraction>0.05 and qivol <= 0.05 and qivol>0.001 and np.isclose(QL[k],0) and np.isclose(QR[k],0) and random.random()< pCloud):
            colIdx[k] = 2
            
        #if we have significant cloud frac and enough stuff for thick cloud
        elif( cloud_fraction > 0.05 and ( ( qivol>0.05) or (qlvol>0.05) ) ):

            #randomly select based on pdf
            if(random.random() < pCloud):
                colIdx[k] = 3
                #Ok, so we have a cloud. Is it game over?
                ip=ix+1
                im=ix-1
                if (ip>ixMax):
                    ip = ixMax
                if(im<0):
                    im = 0
                #Check cloud fraction of neighbors on horizontal if it's persistent >0.9 game over no more returns in column
                if( all(CF[im:ip,k] >0.95)):
                    colIdx[k] = 8
                    #print('game over cloud')
                    break
            elif( (bscat[k]>bmin) and (bscat[k]<bmax) ):
                #print('missed cloud!')
                colIdx[k] = 1
        #if the aerosol is above threshold (too thick) check horizontal for wide area of aerosol outbreak
        if( (bscat[k]>=bmax) and colIdx[k]==0 ):
            ip=ix+1
            im=ix-1
            if (ip>ixMax):
                ip = ixMax
            if(im<0):
                im = 0
            if( all(bbscat[im:ip,k]>=bmax) ):
                colIdx[k] = 5
                print('game over aerosol')
                break
            else:
                colIdx=1
        #if there is a small backscatter check vertical neighbors for low backscatter if any are also low -> too thin 
        if( bscat[k]>1e-8 and bscat[k]<=bmin and colIdx[k]==0):
            ip=k+1
            im=k-1
            if (ip>izMax):
                ip = izMax
            if(im<0):
                im = 0
            if ( all( bscat[im:ip] < bmin) ):
                colIdx[k] = 0
            else:
                colIdx[k] = 4
    return colIdx          
 

if __name__ == '__main__' :
    parser = argparse.ArgumentParser( description = 'read gz text file,select +/-30 deg scan angle, and output to another file using --in directory and --out directory')
    parser.add_argument('--in', help = 'path to ncdiag', required = True, dest = 'path')
    parser.add_argument('--out', help = 'path to ncdiag', required = True, dest = 'outpath')
    arg = parser.parse_args()
    print("Reading {}".format(arg.path))
    main(arg.path,arg.outpath)

