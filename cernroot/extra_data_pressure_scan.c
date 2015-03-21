{
    std::vector<std::vector<double> > xs;
    std::vector<std::vector<double> > ys;

    TFile tf ( "data/pion_peak/data_pressure_scan.root" );
    TCanvas * c;
    tf.GetObject( "c1", c );
    TObject *obj; 
    TIter inext(c->GetListOfPrimitives());
    while ((obj=inext())) {
        cout << "Reading: " << obj->GetName() << endl;
        if (obj->InheritsFrom("TGraph")) {
            TGraph * g = static_cast<TGraph*>( obj );
            cout << "graph: " << obj->GetName() << endl;
            int n = g->GetN();
            double * px = g->GetX();
            double * py = g->GetY();

            std::vector<double> x( px, px + n );
            std::vector<double> y( py, py + n );

            xs.push_back( x );
            ys.push_back( y );
        }
    }  

    std::ofstream ofs( "data/pion_peak/data_pressure_scan.dat" );

    std::ostream& os( ofs );

    int nl = xs.size();
    assert ( nl == ys.size() );

    int np = xs[0].size();

    for ( int p = 0 ; p != np ; ++p )
    {
        os << xs[0][p]*1000  ;
        for ( int l = 0 ; l != nl ; ++ l )
        {
            os << " " << ys[l][p] ;
        }
        os << std::endl;
    }
}
