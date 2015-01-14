#include <functional>
#include <map>
#include <utility>
bool is_ok( TH1D * h )
{
    double rms = h->GetRMS();
    double rmserr = h->GetRMSError();
    return rms < 0.4 && rmserr < 1e-2;
}

void set_log( TCanvas * c, int arg = 1 )
{

    TObject *obj;
    TIter next(c->GetListOfPrimitives());
    while ((obj = next())) {
        if (obj->InheritsFrom(TVirtualPad::Class())) 
        {
            static_cast<TCanvas*>(obj)->SetLogy( arg );
        }
    }
}

int count_pads( TCanvas * c )
{
    int count = 0;
    TObject *obj;
    TIter next(c->GetListOfPrimitives());
    while ((obj = next())) {
        if (obj->InheritsFrom(TVirtualPad::Class())) 
        {
            ++count;
        }
    }
    return count;
}


std::vector<int> pmts =
{
    15,
    22, 23, 24, 25, 26, 27 ,
    32, 33, 34, 35, 36, 37, 38,
    41, 42, 43, 44, 45, 46, 47, 48,
    51, 52, 53, 54, 55, 56, 57, 58, 59,
    61, 62, 63, 64, 65, 66, 67, 68,
    72, 73, 74, 75, 76, 77, 78,
    82, 87
};

void fdo( std::function<void(TH1D*, int)> f, TFile * tf, bool newpmts = false  )
{

    std::vector<int> new_pmts = 
    { 15,
        23, 24, 25, 26,
        33, 34, 35, 36, 37, 
        42, 43, 44, 45, 46, 47, 
        53, 54, 55, 56, 57,
        62, 63, 64, 65, 66, 67,
        73, 74, 75, 76, 77 };

    std::vector<int> all_pmts =
    {
        15,
        22, 23, 24, 25, 26, 27 ,
        32, 33, 34, 35, 36, 37, 38,
        41, 42, 43, 44, 45, 46, 47, 48,
        51, 52, 53, 54, 55, 56, 57, 58, 59,
        61, 62, 63, 64, 65, 66, 67, 68,
        72, 73, 74, 75, 76, 77, 78,
        82, 87
    };

    std::vector<int>& pmts = newpmts ? new_pmts : all_pmts;

    for ( int sector : std::vector<int>{ 1,2,3,4,5,6,7,8} )
    {
        for ( int pmt : pmts )
        {
            TH1D * h ;
            tf->GetObject( Form("PmtTiming/hPMT%d", 100*sector + pmt ), h );

            if ( h == 0 )
            {
                std::cout << "missing " << pmt << std::endl;
                continue;
            }

            f( h, 100*sector + pmt );
        }
    }
}

TCanvas * get_canvas( TCanvas * nc )
{
    if ( nc == 0 )
    {
        nc =  new TCanvas( "c", "c");
        return nc;
    }
    else
    {
        return nc;
    }
}

TH1D * extract( TFile *tf, bool newpms = false,  double maxrms = 0.5, double maxrmerr = 1e-2)
{
    TH1D * hdummy =  0;
    tf->GetObject(  "PmtTiming/hPMT115", hdummy );
    TH1D * hsum = static_cast<TH1D*>( hdummy->Clone("hsum" ));

    auto f = [hsum]( TH1D * h , int pmt )
    {
        if ( is_ok( h) ) { hsum->Add( h ); }
    };

    fdo( f, tf, newpms );
    return hsum;
}

TH1D * rms_plot( TFile* tf, TCanvas * nc = 0  )
{
    TCanvas * c  = get_canvas( nc );
    TH1D * hrms  = new TH1D( "hrms", "PMT TIME RMS", 100, 1, -1 );
    auto f = [hrms]( TH1D * h, int pmt ){ hrms->Fill( h->GetRMS() ); };

    fdo( f, tf );

    c->cd();
    hrms->Draw();
    return hrms;
}

TH1D * fit_centre_plot( TFile* tf, TCanvas * nc = 0 ,
        bool pauseall = false, bool pausethresh = false )
{
    TCanvas * c  = get_canvas( nc );
    TH1D * hcentre  = new TH1D( "hcentre", "PMT TIME CENTRES", 100, 1, -1 );
    std::map<int, double> t0s;
    gStyle->SetOptStat(111111);

    auto f = [hcentre,c, pauseall, pausethresh, &t0s]( TH1D * h, int pmt ) mutable{ 
        if ( is_ok( h ) )
        {
            Double_t c0   = h->GetBinCenter(h->GetMaximumBin());
            Double_t cmin = c0 - 0.3;
            Double_t cmax = c0 + 0.3;

            char name[100];
            TF1 *fit = new TF1 ("fit", "gaus", cmin, cmax);
            h->Fit(fit, "", "", cmin, cmax);
            double fNewT0 = fit->GetParameter(1);

            if ( pauseall || ( pausethresh && fabs(fNewT0) > 0.1 ) )
            {
                h->Draw();
                c->Update();
                std::cout << "Fit is " << fNewT0 << ".\nContinue?" ;
                std::cin.ignore(); //why read something if you need to ignore it? :)
                std::cout << std::endl;
            }
            hcentre->Fill( fNewT0) ;
            t0s.insert( std::make_pair( pmt, fNewT0 ) );
            std::cout << t0s.size() << std::endl;
        }
    };

    std::cout << t0s.size() << std::endl;

    fdo( f, tf );

    std::ofstream fnew_t0s("output/newt0s.dat" );

    for ( auto& pair : t0s )
    {
        fnew_t0s 
            << std::setw(3) << pair.first
            << std::setw(15) << pair.second
            << std::endl;
    }

    c->cd();
    hcentre->Draw();
    return hcentre;
}


TCanvas * plot( TFile* tf, int sector, TCanvas * nc = 0  )
{
    std::cout << "Start" << std::endl;
    TCanvas * c = get_canvas( nc );
    int npads = count_pads( c );
    std::cout << "NPADS" << npads << std::endl;
    if ( count_pads( c ) < 10 )
    {
        c->Divide( 9, 8 );
    }

    auto f = [c, sector]( TH1D* h, int pmt )
    {
        int this_sector = pmt / 100;
        if( this_sector != sector ) return;

        int sector_pmt = pmt - ( 100* this_sector );

        int row = sector_pmt / 10;
        int column = sector_pmt % 10;

        //std::cout << pmt << " " << sector << " " << column << " " << row << std::endl;

        c->cd( (row-1) * 9 + column );

        if (is_ok( h))
        {
            h->SetLineColor( kGreen +2 );
        }
        else
        {
            h->SetLineColor( kRed +2 );
        }

        h->Draw();
    };

    fdo( f, tf );

    return c;
}


