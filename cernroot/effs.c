////////////////////////////////////////////////////////////////////////////

double nm=1;

// Hamamatsu R7400U-03 quantum efficiency
// Parameterized by Evgueni, July 2011

Double_t QE_R7400U_03 (Double_t wavelength) {
    Double_t par[10] =
    {-58.41145755814, 1.450540667766, -0.01561331198442,
        9.545010080831e-05, -3.648461145542e-07, 9.047599515597e-10,
        -1.457151808585e-12, 1.471328774241e-15, -8.46121819724e-19,
        2.11384701372e-22};

    Double_t x = wavelength/nm;
    if (x<180) return 0;
    if (x>660) return 0;
    Double_t qe = 0;
    for (int i=0; i<10; i++) qe += par[i]*TMath::Power(x,i);
    if (qe<0) qe = 0;
    return qe;
}

// Hamamatsu R9880U-110 quantum efficiency
// Parameterized by Angela, July 2011

Double_t QE_R9880U_110 (Double_t wavelength) {
    Double_t par[15] =
    {1489.053765000671, -20.61340505642701, 0.09607362916193821,
        -0.000144918944048782, -1.087924475686453e-07, 3.619104979507752e-10,
        2.742092765095943e-13, -1.067200613381487e-15, 6.333980140159196e-19,
        4.675391577876988, 505.1903283978535, 15.37334879108591,
        -23.08738129086531, 358.7521218115685, 53.63424346389683};

    Double_t x  = wavelength;
    Double_t x1 = (x<650) ? x : 650;
    Double_t qe = 0;
    for (int i=0; i<9; i++) 
    {
        std::cerr << qe << std::endl;
        qe += par[i]*TMath::Power(x1,i);
    }
    std::cerr << qe << std::endl;
    qe += par[9]*TMath::Gaus(x1, par[10], par[11]);
    std::cerr << qe << std::endl;
    qe += par[12]*TMath::Gaus(x1, par[13], par[14]);
    std::cerr << qe << std::endl;
    qe *= 0.01;
    if (x>650) qe *= (1 - (x-650)/(675-650));
    if (qe<0 || x<200) qe = 0;
    return qe;
}

// Hamamatsu R9880U-210 quantum efficiency
// Parameterized by Evgueni, March 2013

Double_t QE_R9880U_210 (Double_t wavelength) {
    double par[9] = 
    {277.3385690654, -5.360192324445, 0.04415739632667,
        -0.0002031054657811, 5.721437395991e-07, -1.012602804374e-09,
        1.100802213492e-12, -6.72600529683e-16, 1.769806940956e-19};

    Double_t x  = wavelength/nm;
    Double_t x1 = (x<680) ? ((x>240) ? x : 240) : 680;
    Double_t qe = 0;
    for (int i=0; i<9; i++) qe += par[i]*TMath::Power(x1,i);
    if (x>680) qe *= (1 - (x-680)/(700-680));
    if (x<240) qe *= (1 - (240-x)/(240-225));
    if (qe<0)  qe = 0;
    return qe;
}

/** Original Cedar PMT quantum efficiency **/

// 1) Lau's parameterization

Double_t QE_EMI_9820_QB_Lau (Double_t wavelength) {
    Double_t wl = wavelength/nm;
    Double_t qe = 0.25 - TMath::Power((wl-400)/500., 2);
    if (qe<0) qe = 0;
    return qe;
}

// 2) Francis's parameterization following the data sheet (March 2013)

Double_t QE_EMI_9820_QB (Double_t wavelength) {

    Double_t wlraw = wavelength/nm;
    Double_t wl = (wlraw < 131.0 ) ? 131.0 : (wlraw > 659.0) ? 659.0 : wlraw;

    Double_t wls[26] =
    {130, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
        340, 360, 380, 400, 420, 440, 460, 480, 500, 520,
        540, 560, 600, 650, 660};

    Double_t qes[26] =
    {0, 0, 20, 21.6, 21.6, 21.4, 21.4, 22, 24, 25, 25.6,
        26, 26, 26.4, 26, 24.8, 23.2, 20.8, 18, 15.2, 12,
        8, 5.6, 2, 0, 0};

    unsigned int i = 0;
    while (true) {
        if (wl == wls[i] || (wl > wls[i] && wl < wls[i+1]) )
        { break; }
        else
        { ++i; }
    }
    return 0.01 * (qes[i] + ( wl - wls[i] )
            / ( wls[i+1] - wls[i] ) * ( qes[i+1] - qes[i] ) );
}


Double_t QE (Double_t wavelength, Int_t PMType) {
    if (PMType==1) return QE_EMI_9820_QB (wavelength);
    if (PMType==2) return QE_R7400U_03   (wavelength);
    if (PMType==3) return QE_R9880U_110  (wavelength);
    if (PMType==4) return QE_R9880U_210  (wavelength);
    return 0;
}

//Neon: parametrization from http://refractiveindex.info // T=273K
Double_t NeRefIndex (Double_t wavelength, Double_t PressureFactor, Double_t TemperatureFactor)
{
    Double_t x            = wavelength;
        Double_t ReducedIndex = 0.012055 * (0.1063*x*x/(184.661*x*x-1) + 1.8290*x*x/(376.840*x*x-1));
    ReducedIndex         *= (PressureFactor/TemperatureFactor);
    Double_t Index        = 1 + ReducedIndex;
    return Index;
}

// Nitrogen: parametrization from http://refractiveindex.info
// T=273K
Double_t N2RefIndex (Double_t wavelength, Double_t PressureFactor, Double_t TemperatureFactor) {
    Double_t x            = wavelength;
    Double_t ReducedIndex = 68.5520E-6 + 32431.57E-6*x*x/(144*x*x-1);
    ReducedIndex         *= (PressureFactor/TemperatureFactor);
    Double_t Index        = 1 + ReducedIndex;
    return Index;
}

void  n2()
{
    std::cerr << N2RefIndex( 0.5, 1.710 / 1.01325, 293.15/273.15 ) << std::endl;
}

void run( std::string output )
{
    std::cerr << QE( 400, 3 ) << std::endl;
    return 0;

    std::streambuf * buf;
    std::ofstream of;

    if ( output != "cout" )
    {
        of.open( output );
        buf = of.rdbuf();
    }
    else
    {
        buf = std::cout.rdbuf();
    }

    std::ostream os ( buf );

    os 
        << std::setw(15) << "wavelength"
        << std::setw(15) << "QE_EMI_9820_QB"
        << std::setw(15) << "QE_R7400U_03"
        << std::setw(15) << "QE_R9880U_110"
        << std::setw(15) << "QE_R9880U_210"
        << "\n";

    for ( int l = 0 ; l != 1000 ; ++l )
    {
        os << std::setw(15) << l;

        for ( int pmt = 1 ; pmt != 5 ; ++pmt )
        {
            os << std::setw(15) << QE( l, pmt);
        }

        os << "\n";
    }
}
