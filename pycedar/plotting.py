import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def blind_colormap( source, name, blind):
    brbgcmp = plt.get_cmap( source )
    brbgdata = brbgcmp._segmentdata.copy()

    for colour in ['red', 'green', 'blue' ]:
        #jiggery
        colorlist = brbgdata[colour]
        clist = [  list( i ) for i in colorlist  ]
        clist[5][1] = clist[4][2]
        clist[5][2] = clist[6][1]
        clist.insert( 6, clist[5] )
        clist[5] = [ 0.5 - 0.5*blind, clist[5][1], 0.7 ]
        clist[6] = [ 0.5 + 0.5*blind, 0.7 , clist[6][2]  ]
        colorlist = [ tuple( i ) for i in clist ]
        brbgdata[colour] = colorlist

    newcm = mpl.colors.LinearSegmentedColormap( name , brbgdata )
    return newcm

def len_mm( x, pos ):
    return '{0:1.1f}'.format( x*1e-3)

def len_mmm( x, pos ):
    return '{0:1.2f}'.format( x*1e-3)

format_mm = FuncFormatter( len_mm )

format_mmm = FuncFormatter( len_mmm )

def plot_dataset( ax, dt ):
    return ax.bar( dt.index, dt.values)
