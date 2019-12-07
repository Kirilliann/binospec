from regions import read_ds9, write_ds9
import matplotlib
matplotlib.use('pdf')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import os
import json
import numpy as np
import glob
from PIL import Image

def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])

def make_thumbnail(ra, dec, outdir, out_fname_img, outfname, region_list, 
                   aux_regs=[], colors=[], pscale=0.27, size=200, imglayer='decals-dr7',
                   compass=False, compass_length = (10,10), pa_out=False):
        fits_files = glob.glob(os.path.join(outdir,'*.fits'))
        jpeg_files = glob.glob(os.path.join(outdir,'*.jpeg'))
        jpeg_command = "wget 'http://legacysurvey.org/viewer/jpeg-cutout/?ra=%.5f&dec=%.6f&layer=%s&pixscale=%.2f&bands=grz&size=%i' -O '%s.jpeg'" %(ra, dec, imglayer, pscale,  size, os.path.join(outdir,out_fname_img))
        if os.path.join(outdir,out_fname_img+'.jpeg') not in jpeg_files:
            os.system(jpeg_command)
        fits_command = "wget 'http://legacysurvey.org/viewer/fits-cutout/?ra=%.5f&dec=%.6f&layer=%s&pixscale=%0.2f&bands=grz&size=%i' -O '%s.fits'" %(ra, dec, imglayer, pscale, size, os.path.join(outdir,out_fname_img))
        if os.path.join(outdir,out_fname_img+'.fits') not in fits_files:
            os.system(fits_command)
        hdr = fits.open(os.path.join(outdir,out_fname_img)+'.fits')[0].header
        img = mpimg.imread(os.path.join(outdir,out_fname_img)+'.jpeg')
        #img = rgb2gray(img)
        px_size = np.sqrt(hdr['CD1_1']**2 + hdr['CD1_2']**2)
        w = WCS(hdr)
	fig = plt.figure(frameon=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.set_size_inches(size/100.0,size/100.0)
        fig.add_axes(ax)
        ax.imshow(img, cmap='gray', origin='upper')
        for regnm in aux_regs:
            regs = read_ds9(regnm)
            for reg in regs:
                pixel_region = reg.to_pixel(w)
                pixel_region.vertices.y = ax.get_ylim()[0] - pixel_region.vertices.y
                pixel_region.plot(ax=ax, color='green')
        for i,reg in enumerate(region_list):
            pixel_region = reg.to_pixel(w)
            pixel_region.vertices.y = ax.get_ylim()[0] - pixel_region.vertices.y
            color = colors[i] if len(colors)>0 else 'red'
            pixel_region.plot(ax=ax, color=color, linewidth=1)
        fig.savefig(os.path.join(outdir,outfname), dpi=100)
        if pa_out:
            im = Image.open(os.path.join(outdir,outfname))
            im = im.resize((1200,1200), Image.ANTIALIAS)
            im = im.rotate(pa_out)
            #im = im.resize((200,200), Image.ANTIALIAS)
            im.save(os.path.join(outdir,outfname))