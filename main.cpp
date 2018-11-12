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
#include "Test.h"
#include <REDFile/File.hh>

using CLHEP::mV;
using CLHEP::ns;
using namespace std;
using namespace RED;

// PROTOTYPES

// The structure contains parameters for auto (random) set vector of photon arrival times
// Parameter "number" is also logic: if number == 0, vector is set manually
struct ArrTimes {
	int number;      // Number of photons (manual set if 0)
	Double_t leftEdge;  //  Left time edge
	Double_t rightEdge; // Right time edge
};

// The enumerator lists possible types of SPE shape
enum {
	kModeNone,  // 0 - No SPE
	kModeF1,    // 1 - TF1 function
	kModeSpline // 2 - TSpline spline
} SPE_Type;

// The union contain pointer to SPE shape in form of TF1 or TSpline
union {
	TF1*      func;
	TSpline3* spline;
} SPE_Shape;
	
// Set vector of photons arrival times from array of double
void SetTseq (ArrTimes AT, vector <Double_t> *Tseq, Double_t *TseqArr, Int_t TseqArrSize);

// Print all photon arrival times
void PrintTseq (const vector <Double_t> *Tseq);

// Create average OutWave for DebugN OutWaves with the same SPE arrival times
void AverageOW (int DebugN, Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay, MakeWave *aex);

int main () {
	cout << "OK, lets go" << endl;
	new TApplication("Canvas", 0, 0);
	Int_t DebugNum  = 1; // Number of OutWaves for averaging by AverageOW

	// Parameters for random set photons arrival times
	ArrTimes AT;
		AT.number = 0;
		AT.leftEdge  = -100*ns;
		AT.rightEdge =  100*ns;

	Double_t TseqArr[] = { -90, -70, -65, -30 , 0 , 5 , 15 , 30, 55, 60}; // Manual values in ns 

	// Set PMT Hamamatsu R11410-20 parameters.
	// Some values are taken from articles C.H.Faham et.al. "Measurements of wavelength-dependent
	// double photoelectron emission..." //JINST 2015 and D.Yu.Akimov et.al. "Performance of
	// Hamamatsu R11410-20 PMTs..." //JINST 2016
	
	// SPE parameters
	Double_t SPE_Width  = 10*ns;    // Width of SPE pulse shape (FWHM for gauss). Doesn't used for TSpline
	Double_t SPE_Xmin   = -50*ns;   // Begin time of domain
	Double_t SPE_Xmax   =  50*ns;   // End time of domain
	Double_t Area_mean  = 10*mV*ns; // SPE pulse area
	Double_t Area_sigma = 2*mV*ns;  // Sigma for spreading SPE area by gauss
	// SPE shape
	SPE_Type            = kModeF1;  // Type of SPE shape:  0 (kModeNone)   - None ,
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
			Double_t splx_bat[] = {-4.4222, -4.1596, -3.8325, -3.5033, -3.2408, -2.9115, -2.5844, -2.2767, -1.9496, -1.8496, -1.8184, -1.5773, -1.2696, -0.94253, -0.65848, -0.63481, -0.59177, -0.57025, -0.50355, -0.41747, -0.32924, 0, 0.32924, 0.41747, 0.50355, 0.57025, 0.59177, 0.63481, 0.65848, 0.94253, 1.2696, 1.5773, 1.8184, 1.8496, 1.9496, 2.2767, 2.5844, 2.9115, 3.2408, 3.5033, 3.8325, 4.1596, 4.4222}; // ns
			Double_t sply_bat[] = {0.1, 1.7722, 2.8957, 3.7818, 4.5096, 5.1505, 5.6331, 6.1948, 6.5983, 5.2296, 3.7027, 2.6583, 2.0096, 2.0887, 2.4131, 3.7027, 4.9052, 5.8704, 6.9227, 5.8704, 4.9052, 4.6679, 4.9052, 5.8704, 6.9227, 5.8704, 4.9052, 3.7027, 2.4131, 2.0887, 2.0096, 2.6583, 3.7027, 5.2296, 6.5983, 6.1948, 5.6331, 5.1505, 4.5096, 3.7818, 2.8957, 1.7722, 0.1}; // mV
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
	// Interaction parameters
	Double_t QE         = 0.3;      // Quantum Efficiency (full)
	Double_t DPE_PC     = 0.225;    // Double Photoelectron Emission probability for PC
	Double_t DPE_1d     = 0.225;    // Double Photoelectron Emission probability for 1dyn
	Double_t QE_1d      = 0.105;    // Quantum Efficiency for 1dyn
	// Geometric parameters
	Double_t Gain_PC_1d = 13;       // Amplification on first gap (PC-1dyn)
	Double_t GF_1d      = 0.1;      // geometric factor (average geom. prob. for a rndm photon from PC to hit 1st dyn)
	Double_t TOFe_PC_1d = 6*ns;     // ToF e- from PC to 1dyn
	Double_t TOFe_mean  = 30*ns;    // ToF e- from PC to anode
	Double_t TOFe_sigma = 3*ns;     // Sigma for spreading ToF PC-anode by gauss
	// Other internal PMT parameters
	Double_t DCR        = 10e6/(1e9*ns);   // Dark count rate, per ns
	Double_t AP_cont    = 0;        // Afterpulsing probability (continuum)
	Double_t AP_peak    = 0;        // Afterpulsing probability (peak)
	
	// Parameters for output waveform
	Double_t Period = 2*ns; // Time between samples of OutWave
	Double_t Gain   = 0.125*mV; // ADC resolution
	Int_t NumSamples    = ceil ((AT.rightEdge + TOFe_mean - AT.leftEdge) / Period); // Samples number in OutWave
	Double_t Delay  = AT.leftEdge; // Delay from "0" of abs.time to "0" sample of OutWave

	// Create Tseq vector of photons arrival times
	vector <Double_t> Tseq;
	SetTseq (AT, &Tseq, TseqArr, sizeof(TseqArr)/sizeof(TseqArr[0]));
	
	// Test of MakePhotons class here
	MakePhotons *Photons = new MakePhotons ();
	Photons->SimulatePhotons(1500, 0.3);
	Tseq = Photons->GetSimPhotonTimes();

	cout << endl << "Program will start" << endl << endl;

	//TApplication *app = 
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
	R11->CalculateParams();

	// Set parameters given above and create OutWave
	MakeWave *MakeWaveObj = new MakeWave; // Create object of MakeWave class
	MakeWaveObj->SetPMT (R11);            // Set used PMT
	MakeWaveObj->SetOutWave (Period, Gain, NumSamples, Delay); // Set OutWave parameters
	MakeWaveObj->SetPhotonTimes (&Tseq);  // Set sequence of photon arrival times
	AverageOW (DebugNum, Period, Gain, NumSamples, Delay, MakeWaveObj);
	R11->Print();
	//PrintTseq(&Tseq);

	Int_t TotalNum = DebugNum*Tseq.size();
	cout << "It were " << TotalNum << " photons that created:" << endl;
	cout << "nothing:\t"     << R11->NumPheHist->GetBinContent(3)*100/TotalNum  << " %" << endl;
	cout << "1 phe in PC:\t" << R11->NumPheHist->GetBinContent(4)*100/TotalNum  << " %" << endl;
	cout << "2 phe in PC:\t" << R11->NumPheHist->GetBinContent(5)*100/TotalNum  << " %" << endl;
	cout << "1 phe in 1d:\t" << R11->NumPheHist->GetBinContent(2)*100/TotalNum  << " %" << endl;
	cout << "2 phe in 1d:\t" << R11->NumPheHist->GetBinContent(1)*100/TotalNum  << " %" << endl;

	// Draw SPE Shape
	R11->DrawShape("user-defined SPE shape;time, ns;Amplitude, MV");
	
	// Draw histogram for number of photoelectrons
	TCanvas *c2 = new TCanvas();
	c2->SetTitle("Number of phe");
	c2->cd();
	c2->SetLogy();
	R11->NumPheHist->Draw();
	R11->NumPheHist->SetTitle("Number of phe resulted (-1 or -2 == 1 or 2 phe from 1st dynode);N(phe);Events");

	// Draw histogram for pulse area
	TCanvas *c3 = new TCanvas();
	c3->SetTitle("Pulse Area");
	c3->cd();
	c3->SetLogy();
	R11->PulseAreaHist->Draw();

	// Draw histogram for amplitude and delay time
	TCanvas *c4 = new TCanvas();
	c4->SetTitle("Distribution of amplitude and time for Photon Pulses");
	c4->Divide(2,1);
	c4->cd(1);
	R11->AmplHist->Draw();
	c4->cd(2);
	R11->TimeHist->Draw();
	
	// Draw histogram for dark pulses time
	TCanvas *c5 = new TCanvas();
	c5->SetTitle("Distribution of amplitude and time for Dark Pulses");
	c5->Divide(2,1);
	c5->cd(1);
	R11->DarkAmplHist->Draw();
	c5->cd(2);
	R11->DarkTimeHist->Draw();
	c5->WaitPrimitive();
	return 0;
}


// FUNCTIONS IMPLEMENTATIONS

// Create average OutWave for DebugN OutWaves with the same SPE arrival times
void AverageOW (int DebugN, Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay, MakeWave* MakeWaveObj) {
	vector <Double_t> MeanOutWave; //Sum of DebugN OutWaves
	MeanOutWave.resize (NumSamples);
	Bool_t ShowAllNum = false; // Report status of each OutWave performing
	                           // If false - only 10 steps will be showed
	// Create OutWave
	if (DebugN < 100) {
		ShowAllNum = true;
	}
	for (int i = 1; i <= DebugN; i++) {
		if (ShowAllNum)
			cout << "Creating " << i << " waveform" << endl;
		else {
			if (!(i % int(floor(DebugN/10)))) {
				cout << "Creating " << i << " waveform" << endl;
			}
		}
		MakeWaveObj->CreateOutWave();
		for (int k = 0; k < NumSamples; k++) {
			MeanOutWave[k] += MakeWaveObj->GetOutWave()[k];
		}
	//cout << i << " OutWave was added to Average Waveform" << endl;
	}
	// Draw OutWave
	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->SetTitle("OutWave");
	TGraph *g2  = new TGraph (NumSamples);
	for (int i = 0; i < NumSamples; i++) {
		g2->SetPoint (i, (Delay + i*Period)/ns, MeanOutWave[i]);
		cout << "point num " << i << "at time= " << (Delay + i*Period)/ns << " ns ";
		cout << "with ampl= " << MeanOutWave[i] << " ADC-units was added" <<endl;
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
			(*Tseq)[i] = AT.leftEdge + rand.Rndm() * (AT.rightEdge - AT.leftEdge);
		}
		sort((*Tseq).begin(), (*Tseq).end());
	}
}

void PrintTseq (const vector <Double_t> *Tseq) {
	cout << "Arrival times:" << endl;
	for (unsigned int i = 0; i < Tseq->size(); i++) {
		cout << ((*Tseq)[i])/ns << " ns" << endl;
	}
}
