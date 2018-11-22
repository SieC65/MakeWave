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

#include "MakeWave.h"
#include "PMT_R11410.hh"
#include "SimPhotons.h"
#include "MakeTest.h"

using CLHEP::mV;
using CLHEP::ns;
using namespace std;
using namespace RED;

int main () {
	cout << "Program start" << endl;
	new TApplication("Canvas", 0, 0);
	
// SET SPE SHAPE

	Double_t SPE_Width  =  10*ns;    // Width of SPE pulse shape (FWHM for gauss). Doesn't used for TSpline
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

	SPE_Type   = kModeF1;  // Type of SPE shape:  0 (kModeNone)   - None ,
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
			Double_t splx[] = {-0.22, 0.05, 0.25, 0.35, 0.61, 0.7, 0.85, 0.89, 0.95}; // ns
			Double_t sply[] = {1    , 2.9 , 5.6 , 7.4 , 9.6 , 8.7, 6.3 , 4.5 , 2   }; // mV
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
	Double_t Area_mean  = 10*mV*ns; // SPE pulse area
	Double_t DCR        = 10e6/(1e9*ns);     // Dark count rate
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
	R11->CalculateParams();
	R11->Print("user-defined");
	R11->Print("calculated");
	R11->Print("probabilities");

// CREATE OUTWAVE
	
	// Set OutWave parameters
	MakeWave *MakeWaveObj = new MakeWave(); // Create object of MakeWave class
	MakeWaveObj->SetPMT (R11);            // Set used PMT
	Double_t Period     = 2*ns;
	Double_t Gain       = 0.125*mV;
	Double_t NumSamples = 150000;
	Double_t Delay      = -150000*ns;
	MakeWaveObj->SetOutWave (Period, Gain, NumSamples, Delay); // Set OutWave parameters

	// Simulate photons
	SimPhotons *Photons = new SimPhotons();
	Int_t NumPhotons = 2000;
	Double_t FracNR     = 0.3;
	Photons->SimulatePhotons(NumPhotons, FracNR);
	cout << "getting simulated photons" << endl;
	vector <Double_t> SimPhotonTimes = Photons->GetSimPhotonTimes();

	cout << "setting simulated photons" << endl;
	MakeWaveObj->SetPhotonTimes (&SimPhotonTimes);  // Set photons arrival times vector

	// Processing OutWave with help of MakeTest class
	cout << "processing average" << endl;
	MakeTest *test = new MakeTest();
	Int_t DebugNum  = 1; // Number of OutWaves for averaging
	test->AverageOW (DebugNum, MakeWaveObj);
	cout << "drawing OW" << endl;
	test->DrawOW(MakeWaveObj);
	cout << "drawing hists" << endl;
	MakeWaveObj->DrawHists();
	// Draw SPE Shape
	test->DrawShape(R11, "user-defined SPE shape;time, ns;Amplitude, MV");
	MakeWaveObj->SaveOutWave("outw.root");

	cout << "well done" << endl;
}