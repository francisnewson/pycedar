import logging
import pandas as pd
import progressbar

def test_aligner( aligner, template_data, test_data ):
    logging.info('Set templates' )
    aligner.set_templates( template_data )

    logging.info('Prepare templates' )
    aligner.prepare_templates()

    logging.info('Prepare datasets' )
    prepared_test_data  = aligner.prepare_data_sets( test_data )

    x_results = []
    y_results = []
    dx_results = []
    dy_results = []

    logging.info('Loop datasets' )
    count = 0
    try:
        with progressbar.ProgressBar( maxval = len(prepared_test_data) ) as pb:
            for pos, dt in aligner.looper(prepared_test_data):
                aligner.compute_alignment( dt )
                best_xy = aligner.best_xy()
                x_results.append( best_xy[0] )
                y_results.append( best_xy[1] )
                dx_results.append( pos[0] - best_xy[0] )
                dy_results.append( pos[1] - best_xy[1] )
                pb.update(count )
                count += 1
    except KeyboardInterrupt:
        pass

    prepared_index = aligner.index( prepared_test_data )[0:count]

    result = pd.DataFrame( {
        'best_x' : pd.Series( x_results, index = prepared_index ),
        'best_y' : pd.Series( y_results, index = prepared_index),
        'dx' : pd.Series( dx_results, index = prepared_index),
        'dy' : pd.Series( dy_results, index = prepared_index) })

    result.reset_index( inplace = True )
    result.columns = ['x', 'y', 'best_x', 'best_y', 'dx', 'dy' ]
    return result
