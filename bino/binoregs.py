import numpy as np
import json
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion, ds9_objects_to_string
from astropy.coordinates import SkyCoord
from plot_regions import make_thumbnail
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument("input_file", help="input config file (json)")
args = parser.parse_args()

hdu = fits.PrimaryHDU()
hdr = hdu.header
json_config = args.input_file #'/home/kirill/Downloads/Coma_dE2b_1_final.json'
data = json.load(open(json_config))
out_reg = open(data["mask_config"]["LABEL"] + ".reg","w+")
theta = -data["mask_config"]["POSA"]/180.0 * np.pi
print data["mask_config"]

for slt in data["sideb"]:
    slt["SIDE"] = "B"

for slt in data["sidea"]:
    slt["SIDE"] = "A"

all_slits = np.concatenate((data["sidea"],data["sideb"]))
rot_mat = np.array([[np.cos(theta), np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
hdr['CD1_1'] = np.cos(theta) * (6.0/3600.0)
hdr['CD1_2'] = np.sin(theta) * (6.0/3600.0)
hdr['CD2_1'] = np.sin(theta) * (6.0/3600.0)
hdr['CD2_2'] = -np.cos(theta) * (6.0/3600.0)
hdr['CRVAL1'] = data["mask_config"]["RA"]
hdr['CRVAL2'] = data["mask_config"]["DEC"]
hdr['CTYPE1'] = 'RA---TAN'
hdr['CTYPE2'] = 'DEC--TAN'
hdr['CRPIX1'] = 0.0
hdr['CRPIX2'] = 0.0
w = WCS(header = hdr)
print w
X_arr = []
Y_arr = []
xdiff = []
ydiff = []

regionslist = []
for obj in all_slits:
    pixc = np.array([[obj['X']],[-obj['Y']]])
    coord0 = np.array([[data["mask_config"]["RA"]],[data["mask_config"]["DEC"]]])
    sky_c = (6.0/3600.0)*np.dot(rot_mat, pixc) 
    sky_c[0] = sky_c[0]  / np.cos(data["mask_config"]["DEC"]/180.0 * np.pi)
    sky_c = sky_c + coord0
    #print obj['RA'] - sky_c[0], obj['DEC'] - sky_c[1]
    mask_coords = np.zeros((4,2))
    polygon = obj['POLYGON']
    for i in range(4):
        mask_coords[i] = w.all_pix2world([[polygon[2*i],polygon[2*i+1]]],1)[0]
    polygon_sky = PolygonSkyRegion(vertices=SkyCoord(mask_coords[:,0], mask_coords[:,1], unit='deg', frame='fk5'))
    regionslist.append(polygon_sky)
    radec = w.all_pix2world([[obj['X'],obj['Y']]],1)
    #print np.sqrt((obj['RA'] - radec[0][0])**2 + (obj['DEC'] - radec[0][1])**2)*3600, obj['RA'], obj['DEC'], mask_coords
    X_arr.append(obj['RA'])
    Y_arr.append(obj['DEC'])
    xdiff.append(obj['RA'] - sky_c[0])
    ydiff.append(obj['DEC'] - sky_c[1])
mask_coords = np.zeros((4,2))
polygon_init = data["mask_config"]["CORNERS"]
polygon = np.array([polygon_init[0],polygon_init[1],polygon_init[0],polygon_init[3],polygon_init[2],polygon_init[3],polygon_init[2],polygon_init[1]])
for i in range(4):
    mask_coords[i] = w.all_pix2world([[polygon[2*i],polygon[2*i+1]]],1)[0]
polygon_sky = PolygonSkyRegion(vertices=SkyCoord(mask_coords[:,0], mask_coords[:,1], unit='deg', frame='fk5'))

regionslist.append(polygon_sky)
mask_coords = np.zeros((4,2))
polygon = np.array([-polygon_init[0],polygon_init[1],-polygon_init[0],polygon_init[3],-polygon_init[2],polygon_init[3],-polygon_init[2],polygon_init[1]])
for i in range(4):
    mask_coords[i] = w.all_pix2world([[polygon[2*i],polygon[2*i+1]]],1)[0]
polygon_sky = PolygonSkyRegion(vertices=SkyCoord(mask_coords[:,0], mask_coords[:,1], unit='deg', frame='fk5'))

regionslist.append(polygon_sky)

out_reg.write(ds9_objects_to_string(regionslist))
out_reg.close()
make_thumbnail(data["mask_config"]["RA"], data["mask_config"]["DEC"], '', data["mask_config"]["LABEL"], data["mask_config"]["LABEL"]+'.png', regionslist, pscale=2.5,size=600, imglayer='dr8')

#plt.plot(Y_arr, ydiff)
#plt.quiver(X_arr, Y_arr, xdiff, ydiff)
#plt.show()
'''
{"RA":195.0902992091987,"WSTART":3786.892896188138,"TARGET":285,"WEND":5301.907412071545,"OBJECT":"J1300
21.67+275354.7","SLIT":1,"pmdec":0,"HEIGHT":0.5009999871253967,"PRIORITY":15,"WIDTH":0.1669999957084656,"epoch":2000,"BBOX":[4.393158781396799,67.55465340194654,0.1669999957084656,0.5009999871
253967],"MAG":16.84461212158203,"OFFSET":0,"Y":67.55465340194654,"X":-51.84244073642547,"DEC":27.89855251926087,"TYPE":"TARGET","pmra":0,"THETA":172.9714584350586,"POLYGON":[-52.10400452855927
,66.11038834641978,-51.87325379328951,67.9819941090375,-51.70625379758104,67.9819941090375,-51.93700453285081,66.11038834641978]}


{"RA":195.0160625,"WSTART":3496.951217677675,"TARGET":583,"WEND":5011.976763057653,"OBJECT":"UDGJ130003.85+280729.3","SLIT":39,"pmdec":0,"HEIGHT":0.5009999871253967,"PRIORITY":1,"WIDTH":0.1669999957084656,"epoch":2000,"BBOX":[-34.99790040594218,-68.49004060226343,0.1669999957084656,0.5009999871253967],"MAG":21.04999923706055,"OFFSET":0,"Y":-68.49004060226343,"X":-91.23349992376444,"DEC":28.124825,"TYPE":"TARGET","pmra":0,"THETA":12.30000019073486,"POLYGON":[-90.82894829786579,-70.72844743142643,-91.37161775382663,-68.23954060870074,-91.20461775811816,-68.23954060870074,-90.66194830215733,-70.72844743142643]}

'''
