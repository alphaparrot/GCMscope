import colormatch
import numpy as np
import matplotlib.pyplot as plt
import sys
import netCDF4 as nc

def orthographic(lon,lat,imap,l0=0,p0=0,ny=36,nx=36,interp='bilinear'):
    xymap = np.zeros((ny,nx,3))
    xymap[:] = 0.0
    coords = np.zeros((ny,nx,2))
    coords2 = np.zeros((ny,nx,2,3))
    coords3 = np.zeros((ny,nx,2))
    rad=0.5*(8*nx/18.0+8*ny/18.0)
    p0 *= np.pi/180.0
    l0 *= np.pi/180.0
    xx = np.arange(nx)-nx/2
    yy = np.arange(ny)-ny/2
    for j in range(0,ny):
        jy = yy[j]
        for i in range(0,nx):
            ix = xx[i]
            if (ix**2+jy**2)<=rad**2:
                rho = np.sqrt(ix**2+jy**2)
                cc = np.arcsin(rho/rad)
                if rho==0:
                    phi = 0.
                    lamb = l0
                else:
                    phi = np.arcsin(np.cos(cc)*np.sin(p0) + jy*np.sin(cc)*np.cos(p0)/rho)
                    lamb = l0 + np.arctan2(ix*np.sin(cc),(rho*np.cos(cc)*np.cos(p0)-jy*np.sin(cc)*np.sin(p0)))
                phi *= 180.0/np.pi
                lamb *= 180.0/np.pi
                if lamb<-180:
                    lamb += 360
                if lamb>180:
                    lamb -= 360
                jlat = np.argmin(abs(phi-lat))
                jlon = np.argmin(abs(lamb-lon))
                if interp=="bilinear":
                    jlat1 = np.where(lat<phi)[0]
                    jlat2 = np.where(lat>=phi)[0]
                    jlon1 = np.where(lon<lamb)[0]
                    jlon2 = np.where(lon>=lamb)[0]
                    if len(jlat1)>0 and len(jlat2)>0 and len(jlon1)>0 and len(jlon2)>0:
                        jlat1 = jlat1[0]
                        jlat2 = jlat2[-1]
                        jlon1 = jlon1[-1]
                        jlon2 = jlon2[0]
                        p1 = lat[jlat1]
                        p2 = lat[jlat2]
                        l1 = lon[jlon1]
                        l2 = lon[jlon2]
                        dl = l2-l1
                        dp = p2-p1
                        if dl==0 and dp>0:
                            fxy1 = imap[jlat1,jlon,:]
                            fxy2 = imap[jlat2,jlon,:]
                            xymap[j,i,:] = (p2-phi)/dp*fxy1 + (phi-p1)/dp*fxy2
                        elif dp==0 and dl>0:
                            fx1y = imap[jlat,jlon1,:]
                            fx2y = imap[jlat,jlon2,:]
                            xymap[j,i,:] = (l2-lamb)/dl*fx1y + (lamb-l1)/dl*fx2y
                        elif dp==0 and dl==0:
                            xymap[j,i,:] = imap[jlat,jlon,:]
                        else:
                            fxy1 = (l2-lamb)/dl*imap[jlat1,jlon1,:] + (lamb-l1)/dl*imap[jlat1,jlon2,:]
                            fxy2 = (l2-lamb)/dl*imap[jlat2,jlon1,:] + (lamb-l1)/dl*imap[jlat2,jlon2,:]
                            xymap[j,i,:] = (p2-phi)/dp*fxy1 + (phi-p1)/dp*fxy2
                    else:
                        jlat1=-200
                        jlat2=-200
                        xymap[j,i,:] = imap[jlat,jlon,:]
                else:
                    xymap[j,i,:] = imap[jlat,jlon,:]
                    
    return xymap

def getphase(phasecurve,nphase):
    ln = phasecurve.variables['lon'][:]
    lt = phasecurve.variables['lat'][:]
    color = phasecurve.variables['colors'][nphase,:,:,:]
    color /= 5.0*np.mean(color)
    proj= orthographic(ln,lt,np.minimum(color,1.0),l0=0.0,nx=200,ny=200)
    return proj

if __name__=="__main__":
    filename = sys.argv[1]
    pc = nc.Dataset(filename,"r")
    tag = filename.split("_phasecurve.nc")[0]
    os.system("mkdir "+tag)
    phases = pc.variables['phase'][:]
    for p in range(0,len(phases)):
        proj = getphase(pc,p)
        f,a=plt.subplots(figsize=(14,12))
        plt.imshow(proj,interpolation='gaussian',origin='lower')
        plt.xticks([])
        plt.yticks([])
        plt.savefig(tag+"/"+tag+"%03d.png",bbox_inches='tight')
        plt.close('all')