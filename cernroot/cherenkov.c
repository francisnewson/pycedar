{
    TFile tfin("data/cherenkov.root");
    TH1D * h;
    tfin.GetObject( "PhotonSpectrum/hPhotonWavelength", h );
    h->Rebin(4);

    unsigned int nBins = h->GetNbinsX();

    for ( unsigned int b = 1 ; b <= nBins ; ++b )
    {
        std::cout
            << std::setw(10) << b
            << std::setw(10) << h->GetBinLowEdge( b )
            << std::setw(10) << h->GetBinContent( b )
            << std::endl;
    }
}
