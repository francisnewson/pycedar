#include "TVirtualFitter.h"
#include "TF1.h"

double gaus(double *x, double *p)
{ return TMath::Gaus(x[0], p[0], p[1]); }

double gausExp( double * x, double * p)
{
    double l = p[0];
    double m = p[1];
    double s = p[2];
    double expPart = TMath::Exp
        ( l/2 * ( 2 * m + l*s*s - 2 * x[0]) );

    double erfPart =
        1 - TMath::Erf( ( m + l*s*s - x[0])
                / (TMath::Sqrt(2) * s ));

    return expPart * erfPart;
}

double reg(double *x, double *p)
{return p[0] * gaus( &(x[0]), &(p[1]));}

double regE(double *x, double *p)
{return p[3] * gaus( &(x[0]), &(p[4]) );}

double regT(double *x, double *p)
{ 
    return p[6] * gausExp ( &(x[0]) , &(p[7]) );
}

double centre_fit(double *x , double *p)
{ 
    if( x[0] > 0.6 && x[0] < 1.7 )
    {
        TF1::RejectPoint();
    }
    return reg(x,p) + regE(x,p) + regT(x,p);
}

double plot_centre_fit(double *x , double *p)
{ 
    return reg(x,p) + regE(x,p) + regT(x,p);
}

double lateE(double *x, double *p)
{
    return p[0] * gaus(x, &p[1]);
}

struct LateFit
{
    TF1 * tfcentre_fit;

    LateFit( TF1 * centre_fit_ )
        :tfcentre_fit( centre_fit_ ){}

    double operator() ( double * x , double * p)
    {
        return plot_centre_fit( x, tfcentre_fit->GetParameters() )
            + lateE( x, p );
    };
};

void fit_2012( TH1D * hdata, bool m_plot = true, bool m_debug = true )
{
    //CENTRE FIT

    //Centre fit function
    TF1 * tf_centre_fit = new TF1( "fCentre", plot_centre_fit, -1, 0.8, 10 );
    tf_centre_fit->SetNpx( 500 );

    tf_centre_fit->SetParameters( 
            1000, 0, 0.2,
            100, -0.5, 0.2,
            100, 0.5, 0.2, 0.05
            );

    double maxdouble = 100000000;
    //Force early and late pulses
    tf_centre_fit->SetParLimits(4 , -10, 0);
    tf_centre_fit->SetParLimits(8, 0, 1);
    tf_centre_fit->SetParLimits(9, 0, 0.1);

    //Ensure positive amplitudes
    tf_centre_fit->SetParLimits(0, 0, maxdouble);
    tf_centre_fit->SetParLimits(3, 0, maxdouble);
    tf_centre_fit->SetParLimits(6, 0, maxdouble);

    tf_centre_fit->SetParNames(
            "AReg", "mReg", "sReg",
            "ARegE", "mRegE", "sRegE",
            "ARegT", "lRegT", "mRegT", "sRegT"
            );

    std::vector<double> x{ 20};
    std::vector<double> p{
        1000, 0, 0.2,
            100, -0.5, 0.2,
            100, 0.5, 0.2, 0.05 };

    //Set up global fitter
    TVirtualFitter *tvf = TVirtualFitter::Fitter(hdata, tf_centre_fit->GetNpar() );
    TVirtualFitter::SetMaxIterations(5000);

    //Do fit
    string fitOpts = "B";
    if (!m_plot) fitOpts += ",N";
    if (!m_debug) fitOpts += ",Q";
    hdata->Fit( tf_centre_fit, fitOpts.c_str() );

    //--------------------------------------------------

    //LATE FIT
#if 0
    auto late_fit = [&tf_centre_fit]( double * x , double * p)
    {
        return centre_fit( x, tf_centre_fit->GetParameters() )
            + lateE( x, p );
    };
#endif

    LateFit * late_fit = new LateFit( tf_centre_fit );


    std::vector<double> xl{ 1.5};
    std::vector<double> pl {100, 1.2, 0.2 };


    TF1 * tf_late_fit = new TF1( "fLate", late_fit, 0.9, 1.6, 3 );

    //Ensure positive amplitudes
    tf_late_fit->SetParLimits(0, 0, maxdouble);

    //Push late means above zero
    tf_late_fit->SetParLimits(1, 1, 2);
    tf_late_fit->SetParLimits(2, 0.1, maxdouble);

    tf_late_fit ->SetParNames(
            "ALateE", "mLateE", "sLateE");

    tf_late_fit->SetParameters(
            100, 1.2, 0.2);

    //Set up global fitter
    TH1D * hdata_late = static_cast<TH1D*>( hdata->Clone( "hdata_late") );
    TVirtualFitter::Fitter(hdata_late, tf_late_fit->GetNpar() );
    TVirtualFitter::SetMaxIterations(5000);

    std::cout << "Doing late fit" << std::endl;
    //Do fit
    hdata_late->Fit( tf_late_fit, fitOpts.c_str());
    std::cout << "Done late fit" << std::endl;

    std::ofstream hfs( "inter/h2012timing.dat" );

    for ( int i = 1 ;  i != hdata->GetNbinsX() ; ++i )
    {
        hfs << hdata->GetBinLowEdge( i ) << " " << hdata->GetBinContent( i ) << std::endl;
    }

    std::ofstream hstats( "inter/h2012stats.dat" );

        for ( int i = 0 ; i != 10 ; ++i )
        {
            hstats << tf_centre_fit->GetParName(i ) <<  " " << tf_centre_fit->GetParameter( i ) << std::endl;
        }

    hstats << std::endl;

    for ( int i = 0 ; i != 3 ; ++i )
    {
        hstats << tf_late_fit->GetParName( i ) << " " <<  tf_late_fit->GetParameter(i) << std::endl;
    }

    tf_late_fit->Draw("same");
}
