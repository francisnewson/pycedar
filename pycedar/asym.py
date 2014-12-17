import math
import pandas as pd

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
