import pandas as pd

def load_templates_csv( f , ignore_pmts = None, index_col = None):
    """Returns a pandas dataframe, extracted from
    a csv file.

    The index is ``r`` ``x`` ``y``, for diaphragm aperture,
    misalignment x and y.

    :param f: csv file containing templates
    :type f: file object
    :rtype: pandas.DataFrame

    """
    if not index_col:
        index_col = ['r', 'x', 'y' ]

    hit_data =  pd.read_csv( f, delim_whitespace = True, index_col = index_col )


    #ensure position names are integers
    hit_data.columns = map( int, hit_data.columns.values )

    if ignore_pmts:
        for pmt in ignore_pmts:
            hit_data[pmt] = 0

    return hit_data

def load_test_sheet( f ):
    """Returns a pandas dataframe, extracted from a csv file
    describing PMT properties

    :param f: csv file containing templates
    :type f: file object
    :rtype: pandas.DataFrame

    """
    pmt_test_sheet = pd.read_csv( f, delim_whitespace=True, index_col = 'serial_number'  )

    #check for duplicates
    assert  pmt_test_sheet.groupby(level=0).filter( lambda x : len(x) > 1 ).empty

    return pmt_test_sheet

def load_pmt_positions( f ):
    """Returns a pandas dataframe, extracted from a csv file
    describing PMT locations

    :param f: csv file containing templates
    :type f: file object
    :rtype: pandas.DataFrame

    """
    pmt_pos_map = pd.read_csv(f , delim_whitespace=True, index_col = 'position' )
    pmt_pos_map.sort_index( inplace = True )
    return pmt_pos_map
