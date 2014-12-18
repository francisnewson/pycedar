import asym
import math
import pandas as pd

def get_octant_corrs( octant_hits ):
     result =   octant_hits.mean() / octant_hits
     result.name = 'corrs'
     return result

def get_octant_corr_stat_err( octant_hits ):
    """Returns statistical errors for octant correction
    factor calculation

    :param octant_hits: pandas DataFrame

    """
    corrs = get_octant_corrs( octant_hits )
    corrs_counts =  pd.concat( [corrs, octant_hits] , axis = 1 )
    corrs_counts.columns = [ 'corrs', 'counts' ]

    N = len( octant_hits )

    return corrs_counts.apply( 
            lambda x : math.sqrt( x['corrs']**2  / x['counts'] 
                - x['corrs'] / ( x['counts']**2 ) / N) ,
            axis = 1 )

    #return ( oc*oc / octant_hits 
             #- oc / (octant_hits * octant_hits ) / len(octant_hits ) )


def full_err( h, s, M , N ):
    return ( ( h * s**2 / M  +1 ) * ( 1 / h - 1 / h / N ) )

def get_octant_corr_full_err( octant_hits, scale_error ):
    """Returns errors for octant correction
    factor calculation assuming a constant fraction error on
    the blue sensitivity data

    :param octant_hits: pandas DataFrame

    """

    corrs = get_octant_corrs( octant_hits )
    corrs_counts =  pd.concat( [corrs, octant_hits] , axis = 1 )
    corrs_counts.columns = [ 'corrs', 'counts' ]
    N = len( octant_hits )
    M = 48

    return corrs_counts.apply( 
            lambda x :  math.sqrt( full_err( x['counts'], scale_error, M, N ) ), 
            axis = 1)
