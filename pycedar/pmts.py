"""About pmts.

pmts.types: a map between PMT types and 
likely start of serial number

"""

pmt_types = {
        'R7400U-03' : { 'HB', 'HC' },
        'R9880U-110' : { 'BAC', 'BAD', 'BDC', 'BHA' }
        }

def get_pmt_type( serial_number):
    """Get PMT type from serial number

    :param serial_number: manufactures serial number
    :type serial_number: string
    :rtype: string

    """
    for pmt_type, serials in pmt_types.items():
        for serial in serials:
            if serial_number.startswith( serial ):
                return pmt_type

    raise KeyError( serial_number )

def get_serial_blue_corr( test_sheet):
    """Get map of serial_number to requred blue correction

    :param test_sheet: PMT Test sheet info
    :type test_sheet: pandas.DataFrame

    """
    pmt_type_groups = test_sheet.groupby( get_pmt_type)['cathode_blue_sens']
    return pmt_type_groups.transform( lambda x: x / x.mean() )

def get_pos_blue_corr( pmt_pos, test_sheet ):
    """Get map of position_id to requred blue correction

    :param pmt_pos: PMT position map
    :type test_sheet: pandas.DataFrame
    :param test_sheet: PMT Test sheet info
    :type test_sheet: pandas.DataFrame

    """
    result = pmt_pos.join( get_serial_blue_corr( test_sheet), on = 'serial_number' )
    return  result.rename( columns = {'cathode_blue_sens' : 'blue_corr' } )
