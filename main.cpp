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
#include <TString.h>

#include "MakeWave.h"
#include "PMT_R11410.hh"

using CLHEP::mV;
using CLHEP::ns;
using namespace std;


// PROTOTYPES

// Print sum of DebugN OutWaves
void AverageOW (int DebugN, const MakeWave::OutWavePar &OW, MakeWave *aex);
// Params of auto (random) set of timeseq - vector of photon arrival times
struct ArrTimes {
	int number;      // Number of arr. times. Manual setting if =0
	Double_t lRange; //  Left range of times
	Double_t rRange; // Right range of times
};
// Set vector of photons arrival times from array of double
void SetTseq (ArrTimes AT, vector <Double_t> *Tseq, Double_t *TseqArr, Int_t TseqArrSize);
// Test function for future
void test (MakeWave* creator, PMT_R11410* pmt);

int main () {
	Int_t DebugNum  = 1; // Number of OutWaves for averaging by AverageOW

	ArrTimes AT;    // Params for random set photons arrival times
		AT.number = 10000000;    // Number of photons (0: manual set)
		AT.lRange = -100*ns; // Begin range of time
		AT.rRange =  100*ns; //   End range of time
	Double_t TseqArr[] = { -20 , 0 , 5.5 , 15 , 50}; // Manual values in ns

	// Type of SPE shape
	enum {
		kModeNone,  // 0 - ???
		kModeF1,    // 1 - TF1
		kModeSpline // 2 - TSpline
	} SPE_Type;

	// SPE shape object
	union {
		TF1*      func;
		TSpline3* spline;
	} SPE_Shape;    

	// Set PMT Hamamatsu R11410-20 parameters.
	// Some values are taken from articles C.H.Faham et.al. "Measurements of wavelength-dependent
	// double photoelectron emission..." //JINST 2015 and D.Yu.Akimov et.al. "Performance of
	// Hamamatsu R11410-20 PMTs..." //JINST 2016
	Double_t SPE_Width  = 10*ns;    // Width of SPE pulse shape (FWHM for gauss)
	Double_t SPE_Xmin   = -50*ns;   // Begin time of domain
	Double_t SPE_Xmax   =  50*ns;   // End time of domain
	SPE_Type            = kModeSpline; // Type of SPE shape: 0 - None , 1 - TF1 , 2 - TSpline
	Double_t QE         = 0.3;      // Quantum Efficiency (full)
	Double_t DPE_PC     = 0.225;    // Double Photoelectron Emission probability for PC
	Double_t DPE_1d     = 0;        // Double Photoelectron Emission probability for 1dyn
	Double_t QE_1d      = 0.105;    // Quantum Efficiency for 1dyn
	Double_t Gain_PC_1d = 13;       // Amplification on first gap (PC-1dyn)
	Double_t GF_1d      = 0.1;      // Average geom. prob. for a rndm photon from PC to hit 1st dyn
	Double_t Area_mean  = 10*mV*ns; // SPE pulse area
	Double_t Area_sigma = 1*mV*ns;  //
	Double_t TOFe_PC_1d = 6*ns;     // ToF e- from PC to 1dyn
	Double_t TOFe_mean  = 30*ns;    // ToF e- from PC to anode
	Double_t TOFe_sigma = 3*ns;     //
	Double_t DCR        = 100;      // Dark count rate, Hz
	Double_t AP_cont    = 0;        //
	Double_t AP_peak    = 0;        // Afterpulsing probability, for continuum and peak
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
			Int_t splpoints = 9;
			Double_t splx[] = {-0.22, 0.05, 0.25, 0.35, 0.61, 0.7, 0.85, 0.89, 0.95}; // ns
			Double_t sply[] = {1    , 2.9 , 5.6 , 7.4 , 9.6 , 8.7, 6.3 , 4.5 , 2   }; // mV
			for (int i=0;i<splpoints;i++) {
				splx[i] = splx[i] * ns;
				sply[i] = sply[i] * mV;
			}
			TGraph* splgraph = new TGraph(splpoints,splx,sply);
			SPE_Shape.spline = new TSpline3("Spline shape",splgraph);
			break;
	}

	MakeWave::OutWavePar OWX; // OutWave parameters
		OWX.Period = 2*ns; // Time between samples of OutWave
		OWX.Gain   = 0.125*mV; // Units of ADC
		OWX.Num    = ceil ((AT.rRange + TOFe_mean - AT.lRange) / OWX.Period); // Samples number in OutWave
		OWX.Delay  = AT.lRange; // Delay from "0" of abs.time to "0" sample of OutWave

	// Create Tseq vector of photons arrival times
	vector <Double_t> Tseq;
	SetTseq (AT, &Tseq, TseqArr, sizeof(TseqArr)/sizeof(Double_t));

	cout << endl << "Program will start" << endl << endl;

	//TApplication *app = 
	new TApplication("Canvas", 0, 0);
	PMT_R11410 *R11 = new PMT_R11410;
	R11->SetDefaults();
	R11->SetParams     (QE, Area_mean, DCR, AP_cont);
	R11->SetDPE_PC     (DPE_PC);
	R11->SetDPE_1d     (DPE_1d);
	R11->SetQE_1d      (QE_1d);
	R11->SetGain_PC_1d (Gain_PC_1d);
	R11->SetGF_1d      (GF_1d);
	R11->SetArea_sigma (Area_sigma);
	R11->SetTOFe_PC_1d (TOFe_PC_1d);
	R11->SetTOFe_mean  (TOFe_mean);
	R11->SetTOFe_sigma (TOFe_sigma);
	R11->SetAP_peak    (AP_peak);
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

	// Set parameters given above and create OutWave
	MakeWave *a = new MakeWave; // Create object of MakeWave class
	a->SetPMT (R11);            // Set used PMT
	a->SetOutWave (OWX);        // Set OutWave parameters
	a->SetTimeSeq (&Tseq);      // Set sequence of photon arrival times
	AverageOW (DebugNum, OWX, a);
	
	R11->DrawShape("user-defined SPE shape;time, ns;Amplitude, arb.un.");
	
	Int_t TotalNum = DebugNum*Tseq.size();
	cout << "It were " << TotalNum << " photons that created:" << endl;
	cout << "nothing:\t"     << a->NPhe->GetBinContent(3)*100/TotalNum  << " %" << endl;
	cout << "1 phe in PC:\t" << a->NPhe->GetBinContent(4)*100/TotalNum  << " %" << endl;
	cout << "2 phe in PC:\t" << a->NPhe->GetBinContent(5)*100/TotalNum  << " %" << endl;
	cout << "1 phe in 1d:\t" << a->NPhe->GetBinContent(2)*100/TotalNum  << " %" << endl;
	cout << "2 phe in 1d:\t" << a->NPhe->GetBinContent(1)*100/TotalNum  << " %" << endl;

	TCanvas *c2 = new TCanvas();
	c2->SetTitle("Number of phe");
	c2->cd();
	a->NPhe->Draw();
	a->NPhe->SetTitle("Number of phe resulted (-1 or -2 == 1 or 2 phe from 1st dynode);N(phe);Events");
	
	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->SetLogy();
	a->PA->Draw();
	c3->WaitPrimitive();
	R11->Print();
	test(a, R11);
	return 0;
}


// FUNCTIONS IMPLEMENTATIONS

// Create average OutWave for DebugN OutWaves with the same SPE arrival times
void AverageOW (int DebugN, const MakeWave::OutWavePar& OW, MakeWave* aex) {
	vector <Double_t> MeanOutWave; //Sum of DebugN OutWaves
	MeanOutWave.resize (OW.Num);
	// Create OutWave
	for (int i = 1; i <= DebugN; i++) {
		cout << "Creating " << i << " waveform" << endl;
		aex->CreateOutWave();
		for (int k = 0; k < OW.Num; k++) {
			MeanOutWave[k] += aex->OutWave[k];
		}
	//cout << i << " OutWave was added to Average Waveform" << endl;
	}
	// Draw OutWave
	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->SetTitle("OutWave");
	TGraph *g2  = new TGraph (OW.Num);
	for (int i = 0; i < OW.Num; i++) {
		g2->SetPoint (i, (OW.Delay + i*OW.Period)/ns, MeanOutWave[i]);
	//	cout << "point num " << i << "at time= " << (OW.Delay + i*OW.Period)/ns << " ns ";
	//	cout << "with ampl= " << MeanOutWave[i] << " ADC-units was added" <<endl;
	}
	g2->SetTitle("OutWave;Time, ns;Amplitude, ADC units");
	g2->Draw();
}

// Set vector of photon arrival times
void SetTseq (ArrTimes AT, vector <Double_t> *Tseq, Double_t *TseqArr, Int_t TseqArrSize) {
	if (!AT.number) {
		// Manual set
		(*Tseq).resize (TseqArrSize);
		for (int i = 0; i < TseqArrSize; i++) {
			(*Tseq)[i] = TseqArr[i]*ns;
		}
	}
	else {
		// Random set
		(*Tseq).resize (AT.number);
		TRandom3 rand;
		rand.SetSeed(0);
		for (int i=0; i<AT.number; i++) {
			(*Tseq)[i] = AT.lRange + rand.Rndm() * (AT.rRange - AT.lRange);
		}
		sort((*Tseq).begin(), (*Tseq).end());
	}
}

// Friend test function
void test (MakeWave* creator, PMT_R11410* pmt) {
		// Get calculated PMT parameters
		/*
		Double_t GetProb_C()          const {return fProb_C;}
		Double_t GetProb_C1()         const {return fProb_C1;}
		Double_t GetProb_C2()         const {return fProb_C2;}
		Double_t GetProb_1d()         const {return fProb_1d;}
		Double_t GetArea_1d()         const {return fArea_1d_mean;}
		Double_t GetArea_1d_Sigma()   const {return fArea_1d_sigma;}
		*/
}