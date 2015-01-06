import math
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import numpy as np
import plotting

def get_octant( position ):
    """Get octant from position ID

    :param position: position ID
    :type f: integer
    :rtype: integer

    """
    return int( position ) // 100

def asym( a, b):
    """Returns normalised assymetry of two numbers
    along with Poisson error

    """
    return ( float( a - b ) / float( a + b ),
            2 * math.sqrt( float(a * b) / float( a + b ) **3 ) )

def oct_asym( octant_hits ):
    """Return UD and LR assymetry for 8 octants

    """
    oh = octant_hits

    u = oh[8] + oh[1]
    d = oh[4] + oh[5]
    l = oh[6] + oh[7]
    r = oh[2] + oh[3]

    ud, ud_err = asym( u, d )
    lr, lr_err = asym( l, r )

    return pd.Series( dict(
        ud = ud, ud_err = ud_err,
        lr = lr, lr_err = lr_err ) )

def octant_corrections( ds ):
    hits = ds.loc[0,0].groupby( get_octant ).sum()
    return hits.mean() / hits


def get_splinefit_data( ds, coord, coord_range ):
    level_coord = 'x' if coord == 'y' else 'y'
    cds = ds.xs( 0, level = level_coord, axis = 0 )
    fit_data = cds[slice(*coord_range)]
    return fit_data

class SplineInfo:
    pass

def get_splines( ds, coord, coord_range, asym, asym_err, s = 10 ):
    fit_data = get_splinefit_data( ds, coord, coord_range )
    asym_spline = UnivariateSpline( fit_data.index.values, fit_data[asym], 
            w = 1 / fit_data[asym_err], s = s )

    coord_points = np.linspace( *coord_range, num = 200 )
    spline_values = asym_spline( coord_points )
    spline_data = pd.DataFrame( { coord : coord_points, asym : spline_values } ).sort( asym )
    inverted_spline = UnivariateSpline( spline_data[asym], spline_data[coord], s = s )

    spi = SplineInfo()
    spi.fit_data = fit_data
    spi.spline_data = spline_data
    spi.spline = asym_spline
    spi.invspline = inverted_spline

    return spi

class AsymAligner:
    def __init__(self, wideset = None, spline_range = 2000, spline_smoothing = 30 ):
        self.corrections = octant_corrections( wideset )
        self.spline_range = spline_range
        self.spline_smoothing = spline_smoothing

    def set_templates( self, templates ):
        self.templates = templates

    def prepare_data_sets( self, data_sets ):
        hits = data_sets.groupby( get_octant, axis = 1 ).sum()
        corr_hits =  hits * self.corrections
        result =  corr_hits.apply( oct_asym, axis = 1 )
        return result

    def looper( self, data_sets):
        return data_sets.iterrows()

    def index( self, data_sets):
        return data_sets.index

    def prepare_templates( self ):
        self.template_group_totals = self.prepare_data_sets( self.templates)
        self.xspi =  get_splines( self.template_group_totals,
                'x', ( -self.spline_range, self.spline_range), 'lr', 'lr_err', self.spline_smoothing )
        self.yspi = get_splines( self.template_group_totals,
                'y', ( -self.spline_range, self.spline_range), 'ud', 'ud_err', self.spline_smoothing )

    def compute_alignment( self, test_data ):
        self.last_x = self.xspi.invspline( test_data['lr'] )
        self.last_y = self.yspi.invspline( test_data['ud'] )

    def best_xy( self):
        return ( self.last_x, self.last_y )

    def plot_xinvspline(self, ax ):
        ax.set_xlabel( 'L\R Asymmetry' )
        ax.yaxis.set_major_formatter(plotting.format_mm)
        ax.set_ylabel( 'x (mm)')
        ax.grid( True )
        ax.autoscale()
        ax.set_xlim( -0.1, 0.1 )

        #print( self.xspi.fit_data )
        #print( self.xspi.spline_data )

        ax.errorbar( self.xspi.fit_data['lr'] , self.xspi.fit_data.index,
                xerr = self.xspi.fit_data['lr_err'], fmt='o', color = 'Blue' )
        ax.plot( self.xspi.spline_data['lr'] , self.xspi.spline_data['x'], ls = '-', lw = 2, color = 'Green' )
        asym_range = np.linspace( -0.1,0.1, 200 )
        ax.plot( asym_range, self.xspi.invspline( asym_range), ls = '--', color = 'Red' )

