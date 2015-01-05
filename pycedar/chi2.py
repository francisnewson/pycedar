import math
from operator import add
import numpy as np
import pandas as pd
import itertools

def group_dict( group_list ):
    forward_dict =  enumerate( group_list ) 
    unpacked_keys = sum([[k]*len(v) for k,v in forward_dict], [])
    unpacked_values = sum( group_list, [] )
    reverse_dict = dict( zip( unpacked_values, unpacked_keys ) )
    return reverse_dict

default_errfun = np.sqrt

def chi2_term( data, data_sqerr, mc, mc_sqerr, ratio ):
    return ( ( data - ratio * mc )**2 
            / ( np.sqrt(data_sqerr) + ratio**2 * np.sqrt(mc_sqerr) ) )

def df_chi2( dt, mc ):
    dtsum = float( dt['hits'].sum() )
    mcsum = float( mc['hits'].sum() )
    ratio = dtsum / mcsum
    chi2_terms = chi2_term( dt['hits'], dt['sqerr'], mc['hits'], mc['sqerr'], ratio )
    return sum( chi2_terms )

def extract_group_totals( data, groups, errfun = default_errfun ):
    data_sqerrs = errfun( data )**2
    data_set = pd.DataFrame( { 'hits' : data , 'sqerr' : data_sqerrs } )
    data_totals = data_set.groupby( groups ).sum()
    return data_totals

class Chi2Aligner:
    def __init__(self,  templates = None, groups = None, errfun = default_errfun):
        self.templates = templates
        self.groups = groups
        self.errfun = errfun

    def extract_group_totals(self, data_set ):
        return extract_group_totals( data_set, self.groups, self.errfun )

    def update_template_cache( self ):
        mydict = {}
        for k in self.templates.index:
            mydict[k] = self.extract_group_totals( self.templates.loc[k] )
        self.template_group_totals = pd.Panel( mydict)
        print( self.template_group_totals )

    def loop_alignment( self, dt ):
        results = []
        #print( 'Start loop' )
        dtsum = float( dt['hits'].sum() )
        dthits = dt['hits'].values
        dtsqerr = dt['sqerr'].values
        for k,mc in self.template_group_totals.iteritems():
            ratio = dtsum / mc['hits'].sum()
            mchits = mc['hits'].values
            mcsqerr = mc['sqerr'].values
            chi2_terms = chi2_term( dthits, dtsqerr, mchits, mcsqerr, ratio )
            results.append( sum( chi2_terms ) )
        #print( 'Stop loop' )
        return  pd.Series( results, index = self.template_group_totals.items)

    def panel_alignment( self, test_data ):
        dt = self.extract_group_totals( test_data )
        return self.template_group_totals.apply( lambda x: df_chi2( dt, x ), axis = 'major' )

    def compute_alignment( self, test_data ):
        self.last_result = self.loop_alignment( test_data ).reset_index()
        self.last_result.columns = [ 'x', 'y', 'chi2' ]
        return self.last_result

#  ____
# / ___|_ __ ___  _   _ _ __  ___
#| |  _| '__/ _ \| | | | '_ \/ __|
#| |_| | | | (_) | |_| | |_) \__ \
# \____|_|  \___/ \__,_| .__/|___/
#                      |_|
#
octant_list = range(1, 9 )

octant_thirds_groups =[
        [ 22, 23, 32, 33, 41, 42, 43, 51, 52, 53, 61, 62, 63, 72, 73, 82 ],
        [ 15, 24, 25, 34, 35, 36, 44, 45, 54, 55, 56, 64, 65, 74, 75, 76 ],
        [ 26, 27, 37, 38, 46, 47, 48, 57, 58, 59, 66, 67, 68, 77, 78, 87 ]
        ]

octant_halves_groups =[ [ 15, 22, 23, 24, 25, 26, 27, 32, 33, 34, 35, 36,
            37, 38, 41, 42, 43, 44, 45, 46, 47, 48, 51, 59 ] ,
        [ 52, 53, 54, 55, 56, 57, 58, 61, 62, 63, 64, 65,
            66, 67, 68, 72, 73, 74, 75, 76, 77, 78, 82, 87 ] ]

def get_chi2_sixths():
    sixth_products = itertools.product( *[octant_thirds_groups, octant_halves_groups]  )
    indi_sixths = [ list( set( t[0] ).intersection( set(t[1]) ) ) for t in sixth_products ]
    sixths = [ [o * 100 + n for n in sixth] for sixth in indi_sixths for o in octant_list ]
    return sixths

def get_chi2_octants():
    return [ [o * 100 + n for halve in octant_halves_groups for n in halve ] for o in octant_list ]
