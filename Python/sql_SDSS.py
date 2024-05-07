from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from numpy import *
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
from copy import copy

def sql_query(ra,dec,z1,z2,r):

    query=f"""
        SELECT
        s.objid,
        sz.ra AS ra,
        sz.dec AS dec,
        pz.z AS photoz,
        pz.zerr AS photozerr,
        sz.z AS specz,
        sz.zerr AS speczerr,
        b.distance AS proj_sep,
        s.modelMag_u AS umag,
        s.modelMagErr_u AS umagerr,
        s.modelMag_g AS gmag,
        s.modelMagErr_g AS gmagerr,
        s.modelMag_r AS rmag,
        s.modelMagerr_r AS rmagerr,
        s.type AS obj_type

        FROM
            PhotoObjAll AS s
        JOIN
            dbo.fGetNearbyObjEq({ra}, {dec}, {r}) AS b ON b.objID = s.objID
        JOIN
            Photoz AS pz ON pz.objid = s.objid
        JOIN
            specObjAll AS sz ON sz.bestobjid = s.objid
        WHERE
            s.type = 3 AND
            sz.z > {z1:.6f} AND
            sz.z < {z2:.6f}
            """
    
    return query

def ra_dec(qso):

    file_systems=open(f'Data/IGM_Danforth_Data/Systems/{qso}_igm-systems.txt','r')

    lines=file_systems.readlines()

    ra=lines[14].split('=')
    ra=ra[1].split('/')[0]
    ra=float(ra)
    
    dec=lines[15].split('=')
    dec=dec[1].split('/')[0]
    dec=float(dec)

    return ra,dec

v_z=lambda z : 3e5*(((1+z)**2-1)/((1+z)**2+1))  # v at z
z_v=lambda v : sqrt((1+((v)/3e5))/(1-((v)/3e5)))-1      # z at v

cosmo=FlatLambdaCDM(H0=69.6,Om0=0.3,Tcmb0=2.725)

qso=['3c263', 'pks0637', 'pks0637', 'pg1424', 'pg0003', 'pg0003', 'pg0003', 'pg1216', 's135712', '1es1553', 'sbs1108', 'pg1222', 'pg1116', 'h1821', 'h1821', 'pg1121', 'pks0405', 'he0056', 'pg1216', '3c263', 'pg1222', 'rxj0439', 'uks0242', 'pg1259', 'pks1302', '3c57', 'p1103', 'phl1811', 'pg0832']

z_abs=array([0.140756, 0.161064, 0.417539, 0.147104, 0.347586, 0.386089, 0.421923, 0.282286, 0.097869, 0.187764, 0.463207, 0.378389, 0.138527, 0.170006, 0.224981, 0.192393, 0.167125, 0.043265, 0.006328, 0.063397, 0.054479, 0.005568, 0.06385, 0.046284, 0.094839, 0.07743, 0.003934, 0.080928, 0.017505])

dA_scale=cosmo.kpc_proper_per_arcmin(z_abs)

dsep=10 #Mpc  transverse distance at absorber redshift to lookfor
arcmin=(dsep*1000)/dA_scale    #  angular separation on plane of sky

v_abs=v_z(z_abs)

z1=z_v(v_abs-1000)
z2=z_v(v_abs+1000)


# for i in range(len(qso)):
# for i in range(2,3):

#     ra,dec=ra_dec(qso[i])
#     query=sql_query(ra,dec,z1[i],z2[i],arcmin[i])

#     data=SDSS.query_sql(query)
#     print(data)
#     print(f'{qso[i]} {z_abs[i]} {len(data)}')
#     data.write(f'{qso[i]}_{z_abs[i]}_SDSS_data.fits',overwrite=True)

ra=238.92935000
dec=11.19010278
r=60

query=f"""
        SELECT
        s.objid,
        sz.ra AS ra,
        sz.dec AS dec,
        pz.z AS photoz,
        pz.zerr AS photozerr,
        sz.z AS specz,
        sz.zerr AS speczerr,
        b.distance AS proj_sep,
        s.modelMag_u AS umag,
        s.modelMagErr_u AS umagerr,
        s.modelMag_g AS gmag,
        s.modelMagErr_g AS gmagerr,
        s.modelMag_r AS rmag,
        s.modelMagerr_r AS rmagerr,
        s.type AS obj_type

        FROM
            PhotoObjAll AS s
        JOIN
            dbo.fGetNearbyObjEq({ra}, {dec}, {r}) AS b ON b.objID = s.objID
        JOIN
            Photoz AS pz ON pz.objid = s.objid
        JOIN
            specObjAll AS sz ON sz.bestobjid = s.objid
        WHERE
            s.type = 3 AND
            sz.z > {z1[9]:.6f} AND
            sz.z < {z2[9]:.6f}
            """
# print(query)
data=SDSS.query_sql(query)
print(data)