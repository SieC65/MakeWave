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

namespace CLHEP {
	static const double_t mV = 1.e-3*volt;
}
using CLHEP::mV;
using CLHEP::ns;
using namespace std;


// FUNCTION PROTOTYPES

// Print sum of DebugN OutWaves
void AverageOW (int DebugN, const MakeWave::OutWavePar &OW, MakeWave *aex);
// Parameters for setting timeseq - vector of photon arrival times
struct ArrTimes {
	int number;			// Number of arr. times. Manual setting if =0
	double_t lrange;	// Left range of times
	double_t rrange;	// Right range of times
};
// Set vector of photons arrival times
void SetTseq (ArrTimes AT, vector <double_t> *Tseq, double_t *TseqArr, Int_t TseqArrSize);


int main() {
	Int_t DebugNum	= 1;	// Number of OutWaves for averaging by AverageOW
	
	ArrTimes AT;	// Params for setting photons arrival times
		AT.number = 1000;		// Number of photons (0 - manual set times)
		AT.lrange = -100*ns;	// Begin range of time (for random set)
		AT.rrange =  100*ns;	//   End range of time (for random set)
	double_t TseqArr[] = { -20 , 0 , 5.5 , 15 , 50};	// Manual values in ns

	// Type of SPE shape
	enum {
		kModeNone,	// 0 - ???
		kModeF1,	// 1 - TF1
		kModeSpline	// 2 - TSpline
    } SPE_Type;
	
	// SPE shape object
	union {
		TF1 *f1;
		TSpline *f2;
	} SPE_Shape;	
		
	// Set PMT Hamamatsu R11410-20 parameters. Given for PMT.
	// Some values are taken from articles C.H.Faham et.al. "Measurements of wavelength-dependent
	// double photoelectron emission..." //JINST 2015 and D.Yu.Akimov et.al. "Performance of
	// Hamamatsu R11410-20 PMTs..." //JINST 2016
		double_t SPE_Width 	= 10*ns;	// Width of SPE pulse shape (FWHM for gauss)
		double_t SPE_Xmin	= -50*ns;	// Begin time of domain
		double_t SPE_Xmax	= 50*ns;	// End time of domain
		SPE_Type			= kModeF1;	// Type of SPE form. 0 - None , 1 - TF1 , 2 - TSpline
		double_t QE			= 0.3;	 	// Quantum Efficiency (full)
		double_t DPE		= 0.225; 	// Double Photoelectron Emission P(2phe)/(P(2phe)+P(1phe))
		double_t QE_1d		= 0.105; 	// QE for 1dyn (only 1phe), prob.
		double_t Gain_PC_1d		= 13; 	// Amplification on first gap (PC-1dyn)
		double_t GF_1d		= 0.1;   	// Average geom. prob. for a rndm photon from PC to hit 1st dyn
		double_t Area_mean	= 10*mV*ns;	// SPE pulse area
		double_t Area_sigma	= 1*mV*ns;	//
		double_t TOFe_PC_1d	= 6*ns; 	// ToF e- from PC to 1dyn
		double_t TOFe_mean	= 30*ns; 	// ToF e- from PC to anode
		double_t TOFe_sigma	= 3*ns;		//
		double_t DCR		= 100; 		// Dark count rate, Hz
		double_t AP_cont	= 0;		//
		double_t AP_peak	= 0;		// Afterpulsing probability, for continuum and peak
		switch (SPE_Type) {
			case 0:
				break;
			case 1:
				SPE_Shape.f1	= new TF1("SPE","gaus(0)",SPE_Xmin,SPE_Xmax);
				SPE_Shape.f1->SetParameter(0, 1);	// Amplitude
				SPE_Shape.f1->SetParameter(1, 0*ns);	// Center
				SPE_Shape.f1->SetParameter(2, (SPE_Width)/(2*sqrt(2*log(2))));	// Sigma
				break;
			case 2:
				//SPE_Shape.f2	= new TSpline();
				break;
		}
		
	MakeWave::OutWavePar OWX; // OutWave Parameters. Given for MakeWave
		OWX.Period 		= 2*ns;		// Time between samples of OutWave
		OWX.Gain 		= 1*mV;//0.125*mV;		// Units of ADC
		OWX.Num 		= ceil((AT.rrange + TOFe_mean - AT.lrange)/OWX.Period);	// Number of samples in OutWave
		OWX.Delay 		= AT.lrange;	// Delay from "0" of abs.time to "0" sample of OutWave
				
	// Create Tseq vector of photons arrival times
	vector <double_t> Tseq;
	SetTseq (AT, &Tseq, TseqArr, sizeof(TseqArr)/sizeof(double_t));
	
	cout << endl << "Program will start" << endl << endl;
	
	TApplication *app = new TApplication("Canvas",0,0);
	PMT_R11410 *R11 = new PMT_R11410;
	R11->SetDefaults();
	R11->SetParameters(QE, DPE, QE_1d, Gain_PC_1d, GF_1d, Area_mean, Area_sigma, TOFe_PC_1d, TOFe_mean, TOFe_sigma, DCR, AP_cont, AP_peak);
	R11->DrawShape("default SPE shape;time, ns;Amplitude, arb.un.");
	if (SPE_Type==1) {
		R11->SetShape(SPE_Shape.f1);
	} else if (SPE_Type==2) {
		R11->SetShape(SPE_Shape.f2);
	}
	
	R11->DrawShape("user-defined SPE shape;time, ns;Amplitude, arb.un.");
	R11->Print();

	// Set parameters given above and create OutWave
	MakeWave *a = new MakeWave;	// Create object of MakeWave class
	a->SetPMT (R11);			// Set used PMT
	a->SetOutWave (OWX);		// Set OutWave parameters
	a->SetTimeSeq (&Tseq);		// Set sequence of photon arrival times
	AverageOW(DebugNum, OWX, a);
	
	cout << "It were " << DebugNum*Tseq.size() << " photons that created:" << endl;
	cout << "nothing:\t" 			<< a->NPhe->GetBinContent(8) << endl;
	cout << "1 phe in PC:\t"	<< a->NPhe->GetBinContent(13) << endl;
	cout << "2 phe in PC:\t"	<< a->NPhe->GetBinContent(18) << endl;
	cout << "1 phe in 1d:\t"	<< a->NPhe->GetBinContent(3) << endl;
	TCanvas *c2 = new TCanvas();
	c2->SetTitle("Number of phe");
	c2->cd();
	a->NPhe->Draw();
	a->NPhe->SetTitle("Number of phe resulted (-1 == 1phe from 1st dynode);N(phe);Events");
	c2->WaitPrimitive();
	return 0;
}


// FUNCTIONS IMPLEMENTATIONS

// Create average OutWave for DebugN OutWaves with the same SPE arrival times
void AverageOW (int DebugN, const MakeWave::OutWavePar &OW, MakeWave *aex) {
	vector <double_t> MeanOutWave;		//Sum of DebugN OutWaves
	MeanOutWave.resize(OW.Num);
	// Create OutWave
	for (int i = 1; i <= DebugN; i++) {
		cout << "Creating " << i << " waveform" << endl;
		aex->CreateOutWave();
		for (int k = 0; k < OW.Num; k++) {
			MeanOutWave[k] += aex->OutWave[k];
		}
	cout << i << " OutWave was added to Average Waveform" << endl;
	}
	// Draw OutWave
	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->SetTitle("OutWave");
	TGraph *g2	= new TGraph(OW.Num);
	for (int i = 0; i < OW.Num; i++) {
		g2->SetPoint(i, (OW.Delay + i*OW.Period)/ns, MeanOutWave[i]);
	//	cout << "point num " << i << "at time= " << (OW.Delay + i*OW.Period)/ns << " ns wit ampl= " << MeanOutWave[i] << " ADC-units was added" <<endl;
	}
	g2->SetTitle("OutWave;Time, ns;Amplitude, ADC units");
	g2->Draw();
}

// Set vector of photon arrival times
void SetTseq (ArrTimes AT, vector <double_t> *Tseq, double_t *TseqArr, Int_t TseqArrSize) {
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
			(*Tseq)[i] = AT.lrange + rand.Rndm()*(AT.rrange-AT.lrange);
		}
		sort((*Tseq).begin(), (*Tseq).end());
	}
}