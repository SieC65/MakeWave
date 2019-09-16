#include <iostream>
#include <algorithm>
#include <vector>

#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TSpline.h>
#include <TH2F.h>
#include <TStyle.h>

#include <REDFile/File.hh>
#include <REDEvent/Event.hh>

#include "MakeWave.h"
#include "PMT_R11410.hh"
#include "SimPhotons.h"
#include "MakeTest.h"

using CLHEP::mV;
using CLHEP::ns;
using namespace std;
using namespace RED;

int main () {
	
// SET SPE SHAPE

	// Only TF1 parameters
	Double_t SPE_Width  =  20*ns;    // Width of SPE pulse shape (FWHM for gauss)
	Double_t SPE_Xmin   = -50*ns;    // Begin time of domain
	Double_t SPE_Xmax   =  50*ns;    // End time of domain

	// The enumerator lists possible types of SPE shape
	enum {
		kModeNone,
		kModeF1,
		kModeSpline
	} SPE_Type;

	// The union contain pointer to SPE shape in form of TF1 or TSpline
	union {
		TF1     *func;
		TSpline *spline;
	} SPE_Shape;

	SPE_Type   = kModeSpline;  // Type of SPE shape:  0 (kModeNone)   - None ,
	                       // 1 (kModeF1) - TF1 , 2 (kModeSpline) - TSpline
	switch (SPE_Type) {
		case kModeNone :
			break;
		case kModeF1 :
			SPE_Shape.func = new TF1("SPE","gaus(0)",SPE_Xmin,SPE_Xmax);
			SPE_Shape.func->SetParameter (0, 1*mV); // Amplitude
			SPE_Shape.func->SetParameter (1, 0*ns); // Center
			SPE_Shape.func->SetParameter (2, (SPE_Width)/(2*sqrt(2*log(2)))); // Sigma
			break;
		case kModeSpline :
			// Taken from LED run
			Double_t splx[] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60}; // ns
			Double_t sply[] = {-0.510066, -0.504701, -3.93614, -10.3283, -15.2794, -17.3102, -17.3881, -13.5405, -9.67392, -5.78841, -1.88383, -1.8855, -1.39651, -0.412552, 0.0830006, -0.398009}; // mV
			Int_t splpoints = 0;
			Int_t splXpoints = sizeof(splx)/sizeof(splx[0]);
			Int_t splYpoints = sizeof(sply)/sizeof(sply[0]);
			if (splXpoints > splYpoints)
				splpoints = splYpoints;
			else
				splpoints = splXpoints;
			for (int i = 0; i < splpoints; i++) {
				splx[i] = splx[i] * ns;
				sply[i] = sply[i] * mV;
			}
			TGraph* splgraph = new TGraph (splpoints,splx,sply);
			SPE_Shape.spline = new TSpline3 ("Spline shape",splgraph);
			break;
	}

// SET PMT

	PMT_R11410 *R11 = new PMT_R11410;
	Double_t QE         = 0.3;      // Quantum Efficiency (full)
	Double_t Area_mean  = 60*mV*ns; // SPE pulse area
	Double_t DCR        = 10e4/(1e9*ns);     // Dark count rate
	Double_t AP_cont    = 0;        // Afterpulsing probability (continuum)
	R11->SetParams     (QE, Area_mean, DCR, AP_cont);
	// Interaction parameters
	R11->SetDPE_PC     (0.225);    // Double Photoelectron Emission probability for PC
	R11->SetDPE_1d     (0.225);    // Double Photoelectron Emission probability for 1dyn
	R11->SetQE_1d      (0.105);    // Quantum Efficiency for 1dyn
	// Geometric parameters
	R11->SetGain_PC_1d (13);       // Amplification on first gap (PC-1dyn)
	R11->SetGF_1d      (0.1);      // geometric factor (average geom. prob. for a rndm photon from PC to hit 1st dyn)
	R11->SetTOFe_PC_1d (6*ns);     // ToF e- from PC to 1dyn
	R11->SetTOFe_mean  (30*ns);    // ToF e- from PC to anode
	R11->SetTOFe_sigma (3*ns);     // Sigma for spreading ToF PC-anode by gauss
	// Other internal PMT parameters
	R11->SetAP_peak    (0);        // Afterpulsing probability (peak)
	R11->SetArea_sigma (2*mV*ns);  // Sigma for spreading SPE area by gauss
	switch (SPE_Type) {
		case kModeNone:
			break;
		case kModeF1:
			R11->SetShape (SPE_Shape.func);
			break;
		case kModeSpline:
			R11->SetShape (SPE_Shape.spline);
			break;
	}

	// Define fitting function for SPE Area
	Double_t fitbeg = 60*mV*ns;
	Double_t fitend = 500*mV*ns;
	char fitfunc[256] = "([0]*ROOT::Math::gaussian_pdf(x,[2],[1])+[3]*ROOT::Math::gaussian_pdf(x,[5],[4])+[6]*ROOT::Math::exponential_pdf(x,[7]))/[8]";
	TF1 *SPEAreaPdf = new TF1 ("pdf for SPE Area",fitfunc,fitbeg,fitend);

	// From Rudik presentation (fit for experimental SPE Area distribution)
	// Gaussian 1 : A*exp(-0.5*((x-x0)/sigma)^2))
	Double_t p0 = 3871;               // A
	Double_t p1 = 134.1*mV*ns;        // x0
	Double_t p2 = 34.55*mV*ns;        // sigma
	// Gaussian 2 : A*exp(-0.5*((x-x0)/sigma)^2))
	Double_t p3 = 194;                // A
	Double_t p4 = 287.3*mV*ns;        // x0
	Double_t p5 = 39.16*mV*ns;        // sigma
	// Exponential : exp(A+kx)
	Double_t p6 = 8.351;              // A
	Double_t p7 = -6.687e-3/(mV*ns);  // k

	SPEAreaPdf->SetParameter(0,p0*sqrt(2*3.1415*p2*p2));
	SPEAreaPdf->SetParameter(1,p1);
	SPEAreaPdf->SetParameter(2,p2);
	SPEAreaPdf->SetParameter(3,p3*sqrt(2*3.1415*p5*p5));
	SPEAreaPdf->SetParameter(4,p4);
	SPEAreaPdf->SetParameter(5,p5);
	SPEAreaPdf->SetParameter(6,exp(p6)/-p7);
	SPEAreaPdf->SetParameter(7,-p7);
	SPEAreaPdf->SetParameter(8,1); // Coeff for normalize function
	SPEAreaPdf->SetParameter(8,SPEAreaPdf->Integral(fitbeg,fitend));
	R11->SetPdfAreaSPE (SPEAreaPdf);

	R11->CalculateParams();
	//R11->Print("user-defined");
	//R11->Print("calculated");
	//R11->Print("probabilities");

// CREATE OUTWAVE
	
	// Set OutWave parameters
	MakeWave *MakeWaveObj = new MakeWave(); // Create object of MakeWave class
	MakeWaveObj->SetPMT (R11);            // Set used PMT
	Double_t Period     = 4*ns;
	Double_t Gain       = 0.125*mV;
	Double_t NumSamples = 75000;
	Double_t Delay      = -75000*ns;
	MakeWaveObj->SetOutWave (Period, Gain, NumSamples, Delay); // Set OutWave parameters

	SimPhotons *Photons = new SimPhotons();
	Photons->SetDefFastFract();
	vector <Double_t> SimPhotonTimes;
	OutputFile *outfile;

	TH2F *h_fracER = new TH2F("fracER","",2001,0,2000,101,0,1.01);
	h_fracER->SetMarkerStyle(7);
	h_fracER->SetMarkerColor(4);
	TH2F *h_fracNR = new TH2F("fracNR","",2001,0,2000,101,0,1.01);
	h_fracNR->SetMarkerStyle(7);
	h_fracNR->SetMarkerColor(2);
	Double_t FracTime = 90*ns;
	Double_t Frac = 0;

	TApplication *app = new TApplication("canvas",0,0);

	char name[128];

	outfile = MakeWaveObj->GetNewFile("ER.root");
	for (Int_t NumPhotons = 100; NumPhotons < 4000; NumPhotons += 1) {
		cout << "ER " << NumPhotons << " photons" << endl;
		SimPhotonTimes = Photons->SimulatePhotons(NumPhotons, "ER");
		MakeWaveObj->SetPhotonTimes (&SimPhotonTimes);
		MakeWaveObj->CreateOutWave();
		MakeWaveObj->AddToFile();
		Frac = MakeWaveObj->GetFrac(FracTime);
		if (Frac)
			h_fracER->Fill(MakeWaveObj->GetNumPE(), Frac);
		//sprintf(name, "ER_wf_%d.root", NumPhotons);
		//MakeWaveObj->SaveOutWave(name);
	}
	MakeWaveObj->CloseFile();

	outfile = MakeWaveObj->GetNewFile("NR.root");
	for (Int_t NumPhotons = 100; NumPhotons < 4000; NumPhotons += 1) {
		cout << "NR " << NumPhotons << " photons" << endl;
		SimPhotonTimes = Photons->SimulatePhotons(NumPhotons, "NR");
		MakeWaveObj->SetPhotonTimes (&SimPhotonTimes);
		MakeWaveObj->CreateOutWave();
		MakeWaveObj->AddToFile();
		Frac = MakeWaveObj->GetFrac(FracTime);
		if (Frac)
			h_fracNR->Fill(MakeWaveObj->GetNumPE(), Frac);
		//sprintf(name, "ER_wf_%d.root", NumPhotons);
		//MakeWaveObj->SaveOutWave(name);
	}
	MakeWaveObj->CloseFile();

	TCanvas *c = new TCanvas("c1","",800,600);
	h_fracER->Draw();
	h_fracNR->Draw("SAME");
	h_fracER->SetXTitle("Photoelectrons number");
	h_fracER->SetYTitle("F90");
	gStyle->SetOptStat(1);
	c->SaveAs("F90.root");
	
	cout << "well done" << endl;
}