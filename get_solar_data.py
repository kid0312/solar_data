#!/usr/bin/env python
# -*- coding: utf-8 -*-

import requests,sys
from io import BytesIO
from scipy.ndimage.interpolation import zoom

import numpy as np
import sunpy.map
import astropy.io.fits as fits

#http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_H_AIA-SDO_FITS_Keyword_Document.pdf
##以下のデータ置き場の対応はまた後日
#https://sdo.gsfc.nasa.gov/assets/img/browse/
def get_aia_image_astropy(t,wavelength,size):
    url = 'http://jsoc2.stanford.edu/data/aia/synoptic/{:04}/{:02}/{:02}/H{:02}00/AIA{:04}{:02}{:02}_{:02}{:02}_{:04}.fits'\
        .format(t.year, t.month, t.day,t.hour, t.year, t.month, t.day, t.hour, t.minute, wavelength)
    
    resp = requests.get(url)
    if resp.status_code == 200:
        contents = BytesIO(resp.content)
        #print(contents)
        hdulist = fits.open(contents)
        hdulist.verify('fix')
        img      = hdulist[1].data
        exptime  = hdulist[1].header['EXPTIME']
        length   = hdulist[1].header['NAXIS1']
        acs_mode = hdulist[1].header['ACS_MODE']
        if (acs_mode=='SCIENCE')and(exptime!=0):
            img = zoom(np.array(img),size/length)
            return img
        
    else:
        return np.zeros((size,size),dtype=float)
    
## c_overはcountの上限値　絶対値がこれ以上であればc_overに置き換える    
def get_hmi_image_astropy(t,size,c_over=None):
    url = 'http://jsoc2.stanford.edu/data/hmi/fits/{:04}/{:02}/{:02}/hmi.M_720s.{:04}{:02}{:02}_{:02}{:02}00_TAI.fits'\
    .format(t.year, t.month, t.day, t.year, t.month, t.day, t.hour, t.minute)
    
    resp = requests.get(url)
    if resp.status_code == 200:
        contents = BytesIO(resp.content)
        #print(contents)
        hdulist = fits.open(contents)
        hdulist.verify('fix')
        img      = hdulist[1].data
        length   = hdulist[1].header['NAXIS1']
        if c_over==None:
            ##nan値処理
            img = np.where(np.isnan(img),0.0,img)
            img = zoom(np.array(img),size/length)
            
        else:
            img = np.where(img>np.abs(c_over),np.abs(c_over),img)
            img = np.where(img<-np.abs(c_over),-np.abs(c_over),img)
            ##nan値処理
            img = np.where(np.isnan(img),0.0,img)
            
        return img
        
    else:
        return np.zeros((size,size),dtype=float)
    

def get_aia_image_sunpy(t,wavelength,size):
    url = 'http://jsoc2.stanford.edu/data/aia/synoptic/{:04}/{:02}/{:02}/H{:02}00/AIA{:04}{:02}{:02}_{:02}{:02}_{:04}.fits'\
        .format(t.year, t.month, t.day,t.hour, t.year, t.month, t.day, t.hour, t.minute, wavelength)
    try:
        imap = sunpy.map.Map(url)
        exptime  = imap.meta['EXPTIME']
        acs_mode = imap.meta['ACS_MODE']
        length   = imap.meta['NAXIS1']
        if (acs_mode=='SCIENCE')and(exptime!=0):
            img = zoom(imap.data,size/length)
        return img
    except:
        return np.zeros((size,size),dtype=float)
    
def get_hmi_image_sunpy(t,size,c_over=None):
    url = 'http://jsoc2.stanford.edu/data/hmi/fits/{:04}/{:02}/{:02}/hmi.M_720s.{:04}{:02}{:02}_{:02}{:02}00_TAI.fits'\
    .format(t.year, t.month, t.day, t.year, t.month, t.day, t.hour, t.minute)

    try:
        imap = sunpy.map.Map(url)
        img  = imap.data
        
        if c_over==None:
            ##nan値処理
            img = np.where(np.isnan(img),0.0,img)
            img = zoom(np.array(img),size/length)
            
        else:
            img = np.where(img>np.abs(c_over),np.abs(c_over),img)
            img = np.where(img<-np.abs(c_over),-np.abs(c_over),img)
            ##nan値処理
            img = np.where(np.isnan(img),0.0,img)
            
        return imap
    
    except:
        return np.zeros((size,size),dtype=float)