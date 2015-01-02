void plot_diaphragm_spectrum( TH2D * h )
{
    h->SetTitle("");
    h->SetContour(100);
    h->SetBit( TH1::kNoTitle, true );
    h->Draw("colz");
    h->GetYaxis()->SetTitleOffset(1.5);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetRangeUser(98.05, 101.95 );
    h->Draw("colz");
}

void diaphragm_spectrum()
{
    gROOT->Reset();
    TFile * tfin = new TFile("data/kaon_mcraytracing.root");
    TH2D * hphotons;
    tfin->GetObject( "CedarMCTester/Photons/RayTracing/Diaphragm", hphotons );

    TCanvas c("c","c", 600, 600 );
    c.cd();
    gStyle->SetOptStat(0);
    gStyle->SetGridColor( kGray);
    c.SetRightMargin(0.15);
    c.SetLeftMargin(0.12);
    c.SetGrid();

    plot_diaphragm_spectrum( hphotons);

    c.Print("output/diaphragm_photons.pdf", "pdf" );

    TH2D * hphotoelectrons;
    tfin->GetObject( "CedarMCTester/Photoelectrons/RayTracing/Diaphragm", hphotoelectrons );
    plot_diaphragm_spectrum( hphotoelectrons);

    c.Print("output/diaphragm_photoelectrons.pdf", "pdf" );

    TCanvas * c2 = new TCanvas( "c2", "c2", 900, 300 );
    c2->cd();
    TFile * tfinpion = new TFile("data/pion_mcraytracing.root");
    TH2D * hpepion;
    tfinpion->GetObject( "CedarMCTester/Photoelectrons/RayTracing/Diaphragm", hpepion );

    gROOT->GetColor( kRed+2)->SetAlpha(0.02);
    hphotoelectrons->SetMarkerColor( kRed+2 );
    hphotoelectrons->Draw();
    hphotoelectrons->GetXaxis()->SetRangeUser(98, 104 );
    gROOT->GetColor( kBlue+2)->SetAlpha(0.02);
    hpepion->SetMarkerColor( kBlue+2 );
    hpepion->Draw("same");


    c2->SetGrid();
    c2->Print("output/diaphragm_pionvskaon.pdf", "pdf" );

    TImage * img = TImage::Create();
    img->FromPad(c2);
    img->WriteImage( "output/diaphragm_pionvskaon.png" );

    TH1D * hkaon = hphotoelectrons->ProjectionX( "hkaon", 0, -1 );
    TH1D * hpion = hpepion->ProjectionX( "hpion", 0, -1 );

    TCanvas * c3 = new TCanvas( "c2", "c2", 900, 500 );
    c3->cd();
    gROOT->GetColor( kRed+2)->SetAlpha(1);
    gROOT->GetColor( kBlue+2)->SetAlpha(1);
    hkaon->SetLineColor( kRed+2);
    hkaon->Scale(1.0/11.0);
    hpion->SetLineColor( kBlue+2);
    hpion->GetXaxis()->SetRangeUser(98,104);
    hpion->SetBit( TH1::kNoTitle, true );
    hpion->Draw();
    hkaon->SetBit( TH1::kNoTitle, true );
    hkaon->Draw("same");
    c3->SetGrid();
    c3->Print("output/diaphragm_pkcomp.pdf", "pdf" );
}
