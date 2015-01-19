TCanvas * plotbeam()
{
    TFile * fin = new TFile( "data/beam_profile.root" );

    fin->ls();

    TH1F * hx = 0;
    TH1F * hy = 0;

    fin->GetObject( "hx", hx );
    fin->GetObject( "hy", hy );

    if (( hx == 0 )|| ( hy == 0 ) )
    {
        std::cout << "HELP!" << std::endl;
        return 0;
    }

    gStyle->SetOptStat(0);

    hx->SetTitle( "" );
    hy->SetTitle( "" );

    TCanvas * c = new TCanvas( "c", "c", 500, 500 );

    hx->Draw();
    hy->Draw("same");

    c->SetGrid( 1, 1 );

    gStyle->SetOptStat(0);
    c->Update();

    //c.Print( "output/beam_profile.pdf", "pdf");
    return c;
}

void beam(){ plotbeam() ; }
