{
    std::vector<std::string> names;
    names.push_back( "ym00200" );
    names.push_back( "ym00400" );
    names.push_back( "ym00500" );
    names.push_back( "ym00600" );
    names.push_back( "ym00700" );
    names.push_back( "ym00900" );
    names.push_back( "yp00000" );
    names.push_back( "yp00200" );
    names.push_back( "yp00400" );
    names.push_back( "yp00500" );
    names.push_back( "yp00600" );
    names.push_back( "yp00700" );
    names.push_back( "yp00900" );

    for ( std::vector<std::string>::iterator iname = names.begin(); 
            iname != names.end() ; ++iname )
    {
        std::string inputfile = "data/photons/beam_ray_yscan/cedar_" + *iname + ".root";

        TFile tf( inputfile.c_str() );

        TH1D * h0;
        TH1D * h1;
        TH1D * h2;
        TH1D * h3;

        tf.GetObject( "CedarPhotons/Photoelectrons/RayTracing/DiaphragmR0", h0 );
        tf.GetObject( "CedarPhotons/Photoelectrons/RayTracing/DiaphragmR1", h1 );
        tf.GetObject( "CedarPhotons/Photoelectrons/RayTracing/DiaphragmR2", h2 );
        tf.GetObject( "CedarPhotons/Photoelectrons/RayTracing/DiaphragmR3", h3 );

        if ( ! ( h0 && h1 && h2 && h3 ) )
        {
            std::cout << "eek!" << std::endl;
            exit(1);
        }

        h0->SetLineColor( kOrange + 2 );
        h0->SetLineWidth( 2 );

        h1->SetLineColor( kRed + 2 );
        h1->SetLineWidth( 2 );

        h2->SetLineColor( kMagenta + 2 );
        h2->SetLineWidth( 2 );

        h3->SetLineColor( kBlue + 2 );
        h3->SetLineWidth( 2 );

        gStyle->SetOptStat(0);

        TCanvas c1( "c", "c", 400, 400 );
        TLegend leg( 0.70, 0.6, 0.88, 0.8 );
        h0->SetTitle( iname->c_str() );
        h0->Draw();
        h0->GetYaxis()->SetRangeUser( 0, 270);
        h1->Draw("same");
        h2->Draw("same");
        h3->Draw("same");

        leg.AddEntry( h0, "LB 1", "L" );
        leg.AddEntry( h1, "LB 2", "L" );
        leg.AddEntry( h2, "LB 3", "L" );
        leg.AddEntry( h3, "LB 4", "L" );

        leg.SetBorderSize( 0 );
        leg.Draw("same");
        
        c1.Print( ("output/photons_diaphragm/beam_yscan/" + *iname + ".pdf").c_str(), "pdf");
    }
}
