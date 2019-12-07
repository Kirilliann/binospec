from regions import read_ds9, write_ds9
from plot_regions import make_thumbnail
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import os
import json
import numpy as np
import glob
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument("input_file", help="input config file (json)")
args = parser.parse_args()


json_config = args.input_file #'/home/kirill/Downloads/Coma_dE2a_1_final.json'
data = json.load(open(json_config))


filename = '' #'/home/kirill/Downloads/Coma_dE1a_1.reg'
allcomaregs = glob.glob('*.reg')
aux_regs = []
for creg in allcomaregs:
    if data["mask_config"]["LABEL"] not in creg:
        aux_regs.append(creg)
    else:
        filename=creg

outdir = data["mask_config"]["LABEL"]
if outdir+'/' not in glob.glob('*/'):
    os.system('mkdir '+outdir)
url_prefix = 'http://gal-03.sai.msu.ru/~kirg/'+data["mask_config"]["LABEL"]+'/'
fits_files = glob.glob('*_input.fits')
jpeg_files = glob.glob('*_input.jpeg')
side_b_data = data["sideb"]
wikifname = data["mask_config"]["LABEL"] + "_wikitab.txt"
outjsonfname = data["mask_config"]["LABEL"] + "_data.json"

#out_f = open(os.path.join(outdir,wikifname),"w+")
for slt in data["sideb"]:
    slt["SIDE"] = "B"

for slt in data["sidea"]:
    slt["SIDE"] = "A"

all_slits = np.concatenate((data["sidea"],data["sideb"]))
for obj in all_slits:
        obj_name = obj["OBJECT"]
        slitn = obj["SLIT"]
        side = obj["SIDE"]
        outfname = "SLIT%i%s_%s.png" % (slitn,side,obj_name)
	outfname = outfname.replace(" ", "_")
        out_fname_img = "SLIT%i%s_%s_input" % (slitn,side,obj_name)
	out_fname_img = out_fname_img.replace(" ", "_")
        ra = obj["RA"]
        dec = obj["DEC"]
        regions = read_ds9(filename)
        make_thumbnail(ra, dec, outdir, out_fname_img, outfname, regions, aux_regs=aux_regs, pscale=0.27, size=200, imglayer='dr8')
        #jpeg_command = "wget 'http://legacysurvey.org/viewer/jpeg-cutout/?ra=%.5f&dec=%.6f&layer=decals-dr7&pixscale=0.27&bands=grz' -O '%s.jpeg'" %(ra, dec, out_fname_img)
        #if out_fname_img+'.jpeg' not in jpeg_files:
        #    os.system(jpeg_command)
        #fits_command = "wget 'http://legacysurvey.org/viewer/fits-cutout/?ra=%.5f&dec=%.6f&layer=decals-dr7&pixscale=0.27&bands=grz' -O '%s.fits'" %(ra, dec, out_fname_img)
        #if out_fname_img+'.fits' not in fits_files:
        #    os.system(fits_command)
        #w = WCS(fits.open(out_fname_img+'.fits')[0].header)
        #img=mpimg.imread(out_fname_img+'.jpeg')
        #fig = plt.figure(frameon=False)
        #ax = plt.Axes(fig, [0., 0., 1., 1.])
        #ax.set_axis_off()
        #fig.set_size_inches(4,4)
        #fig.add_axes(ax)
        #ax.imshow(img, aspect='normal', origin='upper')
        #for regnm in aux_regs:
        #    regs = read_ds9(regnm)
        #    for reg in regs:
        #        pixel_region = reg.to_pixel(w)
        #        pixel_region.vertices.y = ax.get_ylim()[0] - pixel_region.vertices.y
        #        pixel_region.plot(ax=ax, color='green')
        #for reg in regions:
        #    pixel_region = reg.to_pixel(w)
        #    pixel_region.vertices.y = ax.get_ylim()[0] - pixel_region.vertices.y
        #    pixel_region.plot(ax=ax, color='red')
        obj["IMGFILE"] = outfname
	out_f = open(os.path.join(outdir,wikifname),"a+")
        out_f.write("|![IMAGE](%s)|%s|%.4f %.4f|%s%s|%.1f|%.1f-%.1f|%i|%.1f\n" % (url_prefix+outfname, obj_name, float(ra), float(dec), 
                                                                             str(slitn), side, float(obj["THETA"]), obj["WSTART"], obj["WEND"],
                                                                             obj["PRIORITY"], obj["MAG"]))
        out_f.close()
	#fig.savefig(outdir+outfname)

fjson = open(os.path.join(outdir,outjsonfname), 'w+')
json.dump(data, fjson)
