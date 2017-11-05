#include <iostream>
#include <iomanip>

#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>

#include "MakeWave.h"

using namespace std;
namespace CLHEP {
	static const double mV = 1.e-3*volt;
}
using CLHEP::ns;
using CLHEP::mV;

// Simple constructor
MakeWave::MakeWave () : fisrnd(true) {
	cout << "MakeWave object was created" << endl;
}

// Set PMT
void MakeWave::SetPMT (PMT *PMTX) {
	fPMT = PMTX;
	cout << "PMT was set" << endl;
}

// Set OutWave parameters
void MakeWave::SetOutWave (const OutWavePar &OWX) {
	fOWX = OWX;
	cout << "OutWave parameters were set" << endl;
}

// Set sequence of photon times
void MakeWave::SetTimeSeq (vector <double> *Tseq) {
	ftimeseq = Tseq;	// Now addresses of ftimeseq and Tseq are equal (pointer Tseq was created in main.cpp)
	cout << "Sequence of SPE's arrival times was set" << endl;
}

// Creating OutWave
void MakeWave::CreateOutWave () {
	cout << "Creating OutWave..." << endl;
	
	Double_t SampleTime = 0;		// Time of sample from "0" of OutWave
	Double_t PulseTime 	= 0;		// Time from "0" of OutWave to "0" of SPE shape
	Double_t PulseAmpl 	= 0;		// Amplitude of SPE shape
	Int_t StartSample 	= 0;		// First sample in SPE domain
	Int_t FinishSample 	= 0;		// Last sample in SPE domain
	Int_t InteractType 	= 0;		// Type of interaction :
									// "0" - no int, "1" - 1phe in PC, "2" - 2phe in PC, "-1" - 1phe in 1d
	
	OutWave.clear();				//  Clear vector OutWave
	OutWave.resize (fOWX.Num, 0);	// Resize vector OutWave
	
	// Some histograms
	// NPhe - number of photons created 0, 1 and 2 phe
	NPhe = new TH1F("NPhe","NumberOfPhe",20,-1.5,2.5);
	// PA - area under pulse SPE (DPE)
	Double_t PA_max = fPMT->GetArea() * fPMT->GetAmplitude() /(ns*mV);	// Maximum pulse area
	PA = new TH1F("PA","PulseArea[mV*ns]",1000,0,PA_max);
	
	// Create PulseArray vector
	RED::PMT_R11410::PulseArray *Pulse = new RED::PMT_R11410::PulseArray;
	
	// Go along all arrival times and fill PulseArray vector
	cout << "Obtaining pulses vector..." << endl;
	for (unsigned int i = 0; i < ftimeseq->size(); i++) {
		InteractType = fPMT->OnePhoton((*ftimeseq)[i], *Pulse);
		NPhe->Fill(InteractType);
	}
	
	// Go along all pulses in PulseArray vector and create OutWave
	cout << "Adding pulses to OutWave..." << endl;
	for (unsigned int i = 0; i < Pulse->size(); i++) {
		// Calculate left and right samples including SPE
		StartSample		=  ceil( (((*Pulse)[i]).fTime - fOWX.Delay + fPMT->GetXmin()) / fOWX.Period );
		FinishSample	= floor( (((*Pulse)[i]).fTime - fOWX.Delay + fPMT->GetXmax()) / fOWX.Period );
		// Get delay time from "0" of OutWave to "0" of SPE shape
		PulseTime		= ((*Pulse)[i]).fTime - fOWX.Delay;
		// Get amplitude of SPE shape
		PulseAmpl		= ((*Pulse)[i]).fAmplitude;
		// Fill pulse area histogram
		PA->Fill( PulseAmpl/mV * fPMT->GetArea()/ns );
		// Limit edges
		if (StartSample < 0)
			StartSample = 0;
		if (FinishSample > fOWX.Num - 1)
			FinishSample = fOWX.Num - 1;
		// Add SPE to OutWave
		for (int s = StartSample; s <= FinishSample; s++) {
			SampleTime = s*fOWX.Period;
			OutWave[s] += PulseAmpl/fOWX.Gain * fPMT->Eval(SampleTime - PulseTime);
		}
	}
	cout << "OutWave was created" << endl;
}

// Print OutWave
void MakeWave::PrintOutWave() {
	cout << "Printing OutWave..." << endl;
	for (int i = 0; i < fOWX.Num; i++) {
		cout << OutWave[i] << "\t(t = " << (fOWX.Delay + i*fOWX.Period)/ns << " ns)" << endl;
	}
	cout << "OutWave was printed" << endl;
}

// Draw OutWave
void MakeWave::Draw () {
	cout << "Drawing OutWave..." << endl;
	TGraph *g1 = new TGraph(fOWX.Num);
	for (int i = 0; i < fOWX.Num; i++) {
		g1->SetPoint(i, (fOWX.Delay + i*fOWX.Period)/ns, OutWave[i]);
	}
	g1->SetTitle("Waveform");
	g1->Draw();
	cout << "OutWave was drawn" << endl;
}