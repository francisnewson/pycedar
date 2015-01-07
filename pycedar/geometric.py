import csv
import pyrr
import math
import pandas as pd
from scipy import optimize
import numpy as np
from matplotlib import _cntr as cntr
import shapely.geometry

def project( pmt ):
    #shift origin
    real_3d_pos = pyrr.Vector3( [pmt['x'], pmt['y'],
        pmt['z']  - 69979.0] )

    #rotate into standard position
    standard_rotation = pyrr.Quaternion.from_z_rotation(
            math.pi / 8.0 + (pmt['sector']-1) * math.pi / 4.0 )

    sp = pyrr.Matrix33( standard_rotation ) * real_3d_pos

    #rearrange coords
    mixed_up_pos = pyrr.Vector3( [sp.x, sp.z +352, sp.y] )

    #repeat rotation
    final_pos = pyrr.Matrix33( standard_rotation ) * mixed_up_pos

    return pd.Series( dict( x = final_pos.x,  y = final_pos.y ) )

def make_geom_chi2_term( x, y, r ):
    def chi2_term( ds):
        return (ds['hits'].values * np.square( 
                r - np.hypot( ( ds['x'].values - x), ( ds['y'].values - y ) ) )).sum()

    return chi2_term

def geom_chi2( xyr, hitmap):
    my_geom_chi2_term = make_geom_chi2_term( *xyr )
    chi2  = my_geom_chi2_term( hitmap )
    return chi2

def circle_fit( hits, projmap ):
    hit_map = projmap.copy()
    hit_map['hits'] = hits
    res = optimize.minimize(  geom_chi2, (0,0, 300), (hit_map,), method = 'Powell' )

    return pd.Series( dict( fitx = res.x[0] , fity = res.x[1], fitr = res.x[2] ) )

class GeomAligner:
    def __init__(self):
        self.pos_map = pd.read_csv( 'params/pmt_3d_pos.dat', delim_whitespace = True,
                names = ['pmt', 'x', 'y', 'z' ], index_col = 'pmt'  )

        self.pos_map['sector']  = self.pos_map.index / 100
        self.proj_map  = self.pos_map.apply( project, axis = 1 )
        self.templates = pd.DataFrame()
        self.template_fits = pd.DataFrame()

    def set_templates( self, templates ):
        self.templates = templates

    def set_template_fits( self, template_fits ):
        self.template_fits = template_fits

    def prepare_data_sets( self, data_sets ):
        result = data_sets.apply( lambda x : circle_fit(x, self.proj_map ), axis = 1 )
        return result

    def looper( self, data_sets):
        return data_sets.iterrows()

    def index( self, data_sets):
        return data_sets.index

    def get_contours( self, df ):
        X = df.columns.values
        Y = df.index.values
        Z = df.values
        x, y = np.meshgrid( X, Y )
        return cntr.Cntr( x, y, Z )

    def prepare_templates( self ):
        if not self.templates.empty and self.template_fits.empty:
            self.template_fits = self.prepare_data_sets( self.templates ).reset_index()
            self.template_fits.to_pickle( 'cache/geo.pck')

        if self.template_fits.empty:
            raise Exception( 'No templates for GeomAligner' )

        fitx = self.template_fits.pivot( 'x', 'y', 'fitx' )
        fity = self.template_fits.pivot( 'x', 'y', 'fity' )
        self.xcont = self.get_contours( fitx.interpolate().bfill().ffill() )
        self.ycont = self.get_contours( fity.interpolate().bfill().ffill() )

    def paths( self, lines ):
        nseg = len(lines)//2
        segs, codes=  lines[:nseg],lines[nseg:]
        seglists = map( lambda x : x.tolist(), segs )
        codelists = map( lambda x : x.tolist(), codes )
        paths =  zip( seglists, codelists )
        return paths

    def compute_alignment( self, test_data ):
        xlines = self.xcont.trace( test_data.fitx )
        ylines = self.ycont.trace( test_data.fity )

        xpaths = self.paths( xlines )
        ypaths = self.paths( ylines )

        intersections = []

        for xpath in xpaths:
            xpoints = xpath[0]
            xcodes = xpath[1]
            xline = shapely.geometry.LineString( xpoints )
            for ypath in ypaths:
                ypoints = ypath[0]
                ycodes = ypath[1]
                yline = shapely.geometry.LineString( ypoints )

                inter = xline.intersection( yline )
                if type( inter ) == shapely.geometry.collection.GeometryCollection:
                    for p in inter:
                        intersections.append( p )
                elif type( inter ) == shapely.geometry.multipoint.MultiPoint:
                    for p in inter:
                        intersections.append( p )
                else:
                    intersections.append( inter )

        self.best_xys = [ (p.x, p.y) for p in intersections ]
        return [ [xpath[0] for xpath in xpaths], [ypath[0] for ypath in ypaths]]


    def best_xy( self ):
        return list(reversed(self.best_xys[0]))
