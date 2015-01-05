import csv

class CircleFitter:
    def __init__(self):
        self.pos_map = {}
        with open('params/pmt_3d_pos.dat') as f:
            reader = csv.DictReader( f, fieldnames = ['pos', 'x', 'y', 'z'] )
            for pmt in reader:
                self.pos_map[pmt['pos']] = self.project( pmt )

    def project( self, pmt ):
        return 
