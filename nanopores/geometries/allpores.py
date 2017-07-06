# (c) 2017 Gregor Mitscha-Baude
"""common interface for all pore geometries

should behave as follows:
    -) defines default parameters for the geometry
    -) provides function that takes geoname, geo params (including dim) plus
       small additional number of non-geo (solver) params
       (h, reconstruct, subs), and returns geo.
    -) optionally, provide functions that returns underlying geometrical object before
       building
    -) allows easy, modular addition of new geometries"""
import numpy as np
from nanopores.tools.utilities import Params, any_params
import nanopores.tools.polygons as polygons
from nanopores.geometries.cylpore import MultiPore

def get_geo(geoname, **params):
    geoclass = geometries[geoname]()
    return geoclass.get_geo(**params)
    
class BasePore(object):
    default = {}
    
    def get_geo(self, h=1., reconstruct=False, **params):
        self.params = Params(self.default, **params)        
        return self.build(h=h, reconstruct=reconstruct)
    
class WeiPore(BasePore):
    
    default = dict(
        R = 100.,
        R0 = 60.,
        H0 = 70.,
        H = 150.,
        x0 = [0, 0, 46],
        rMolecule = 2.1,
        dim = 3,
        no_membrane = True,
        r0 = 13, # pore radius
        angle = 40, # aperture angle in degrees
        lcCenter = 0.3,
        lcMolecule = 0.1,
        h = 10.,
        subs = None,
        reconstruct = False,
        poreregion = True,
    )
    
    def build(self, h=1., reconstruct=False):
        params = self.params
        # SiN membrane thickness (in vertical direction)
        lsin = 50.
        # Au membrane thickness (in vertical direction)
        lau = 40.
        # Au thickness in radial direction
        rlau = 10.
        # SAM layer thickness (in vertical direction)
        lsam = 3
        
        l0 = lau + lsin + lsam
        angle2 = params.angle/2. * np.pi/180.
        tan = np.tan(angle2)
        sin = np.sin(angle2)
        cos = np.cos(angle2)
        
        l = l0/2.
        r0 = params.r0
        r1 = r0 + l0*tan
        rsam = r0 + lsam/cos
        rsin = r0 + lsam/cos + rlau
        R = params.R
        
        sam = [[r0, -l], [r1, l], [R, l], [R, l - lsam],
               [rsam - tan*(lsam - l0), l - lsam], [rsam, -l]]
        au = [sam[5], sam[4], sam[3], [R, -l + lsin], [rsin + tan*lsin, -l + lsin],
              [rsin, -l]]
        sin = [au[5], au[4], au[3], [R, -l]]
        
        p = MultiPore(**params)
        p.add_polygons(sam=sam, au=au, sin=sin)
        
        receptor = polygons.Ball([30.,0.,30.], 7., lc=0.1)
        p.add_balls(receptor=receptor)
        geo = p.build(h, params.subs, reconstruct)
        return geo
    

geometries = dict(
    wei = WeiPore,
    #pughpore = PughPore(),
    #alphahem = AlphaHem,
)

if __name__ == "__main__":
    params = any_params(geoname="wei", h=10.)
    geo = get_geo(**params)
    geo.plot_subdomains()
    geo.plot_boundaries(interactive=True)