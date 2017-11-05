#include <iostream>

#include <TGraph.h>
#include "SystemOfUnits.h"

#include "MakeWave.h"

using std::cout;
using std::endl;
using std::vector;
using CLHEP::ns;

//using CLHEP::mV;

// Simple constructor
MakeWave::MakeWave () {
	PhotonPulse = 0;
	DarkPulse   = 0;
	fPMT        = 0;
	//cout << "MakeWave object was created" << endl;
}

// Set PMT
void MakeWave::SetPMT (RED::PMT* PMTX) {
	fPMT = PMTX;
	//cout << "PMT was set" << endl;
}

// Set OutWave parameters
void MakeWave::SetOutWave (const OutWavePar& OWX) {
	fOWX = OWX;
	//cout << "OutWave parameters were set" << endl;
}

// Set sequence of photon times
void MakeWave::SetTimeSeq (vector <double>* Tseq) {
	ftimeseq = Tseq;
	//cout << "Sequence of SPE's arrival times was set" << endl;
}

// Creating OutWave
void MakeWave::CreateOutWave () {
	//cout << "Creating OutWave..." << endl;

	OutWave.clear();              //  Clear vector OutWave
	OutWave.resize (fOWX.Num, 0); // Resize vector OutWave
	
	// Create PulseArray vector
	if (!PhotonPulse) {
		PhotonPulse = new RED::PMT::PulseArray;
	}
	else PhotonPulse->clear();

	// Create DarkPulseArray vector
	if (!DarkPulse) {
		DarkPulse = new RED::PMT::PulseArray;
	}
	else DarkPulse->clear();

	/// ADD SPE FROM PHOTONS
	//cout << "adding photons" << endl;
	// Generate vector Pulse of PulseArray type
	for (unsigned int i = 0; i < ftimeseq->size(); i++) {
		fPMT->OnePhoton (&((*ftimeseq)[i]), *PhotonPulse, false);
	}
	AddPulseArray (PhotonPulse);

	/// ADD DARK COUNTS
	//cout << "adding dark" << endl;
	// Generate vector DarkPulse of PulseArray type
	fPMT->GenDCR (fOWX.Delay, fOWX.Delay + fOWX.Num * fOWX.Period, *DarkPulse);
	AddPulseArray (DarkPulse);
	//cout << "OutWave was created" << endl;
}

// Add PulseArray vector to OutWave
void MakeWave::AddPulseArray (RED::PMT::PulseArray *Pulses) {
	Double_t SampleTime   = 0; // Time of sample from "0" of OutWave
	Double_t PulseTime    = 0; // Time from "0" of OutWave to "0" of SPE shape
	Double_t PulseAmpl    = 0; // Amplitude of SPE shape
	Int_t StartSample     = 0; // First sample in SPE domain
	Int_t FinishSample    = 0; // Last sample in SPE domain
	
	// Go along all pulses in PulseArray vector DarkPulse and add them to OutWave
	for (unsigned int i = 0; i < Pulses->size(); i++) {
		// Calculate left and right samples including SPE
		StartSample     =  ceil( (((*Pulses)[i]).fTime - fOWX.Delay + fPMT->GetXmin()) / fOWX.Period );
		FinishSample    = floor( (((*Pulses)[i]).fTime - fOWX.Delay + fPMT->GetXmax()) / fOWX.Period );
		// Get delay time from "0" of OutWave to "0" of SPE shape
		PulseTime       = ((*Pulses)[i]).fTime - fOWX.Delay;
		// Get amplitude of SPE shape
		PulseAmpl       = ((*Pulses)[i]).fAmpl;
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