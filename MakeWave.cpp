#include <iostream>

#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
//#include <TString.h>
#include "SystemOfUnits.h"

#include "MakeWave.h"

using std::cout;
using std::endl;
using std::vector;
using CLHEP::ns;

//using CLHEP::mV;

// Simple constructor
MakeWave::MakeWave () {
	fPhotoElectrons   = 0;
	fDarkElectrons    = 0;
	fPMT              = 0;
	fDarkTimeHist     = 0;
	fDarkAmplHist     = 0;
	fNumPheHist       = 0;
	fPulseAreaHist    = 0;
	//cout << "MakeWave object was created" << endl;
}

// Set PMT
void MakeWave::SetPMT (RED::PMT* PMT) {
	fPMT = PMT;
	//cout << "PMT was set" << endl;
}

// Set OutWave parameters
void MakeWave::SetOutWave (Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay) {
	fPeriod     = Period;
	fGain       = Gain;
	fNumSamples = NumSamples;
	fDelay      = Delay;
	//cout << "OutWave parameters were set" << endl;
}

// Set sequence of photon times
void MakeWave::SetPhotonTimes (vector <double>* PhotonTimes) {
	fPhotonTimes = PhotonTimes;
	//cout << "Sequence of SPE's arrival times was set" << endl;
}

// Set some default parameters
void MakeWave::SetDefaults(){
	fPMT = new RED::PMT_R11410();
	fPMT->SetDefaults();
	fPeriod     = 2*ns;
	fGain       = 0.125*mV;
	fNumSamples = 150000;
	fDelay      = -150000*ns;
}

// Creating OutWave
void MakeWave::CreateOutWave () {
	//cout << "Creating OutWave..." << endl;

	fOutWave.clear();              //  Clear vector OutWave
	fOutWave.resize (fNumSamples, 0); // Resize vector OutWave
	Char_t NumPhe = 0;

// ADD SPE FROM PHOTONS

	// Check if PulseArray vector fPhotoElectrons exists
	if (!fPhotoElectrons) {
		fPhotoElectrons = new RED::PMT::PulseArray;
	}
	else fPhotoElectrons->clear();

	// Histogram for number of photons created 0, 1 and 2 phe
	// -2: 2phe from 1d; -1: 1phe from 1d; 0: no interaction
	// +1: 1phe from PC; +2: 2phe from PC
	if (!fNumPheHist) {
		fNumPheHist   = new TH1F ("fNumPheHist","NumberOfPhe",5,-2.5,2.5);
	}

	// Generate fPhotoElectrons, add them to OutWave and fill the phe number hist
	for (unsigned int i = 0; i < fPhotonTimes->size(); i++) {
		NumPhe = fPMT->OnePhoton (&(fPhotonTimes->at(i)), *fPhotoElectrons, false);
		fNumPheHist->Fill(NumPhe);
	}
	AddPulseArray (fPhotoElectrons);

// ADD DARK COUNTS

	// Check if PulseArray vector fDarkElectrons exists
	if (!fDarkElectrons) {
		fDarkElectrons = new RED::PMT::PulseArray;
	}
	else fDarkElectrons->clear();

	// Generate dark electrons and add them to OutWave
	fPMT->GenDCR (fDelay - (fPMT->GetXmax() - fPMT->GetXmin()), fDelay + fNumSamples * fPeriod, *fDarkElectrons);
	AddPulseArray (fDarkElectrons);

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
		StartSample     =  ceil( (((*Pulses)[i]).fTime - fDelay + fPMT->GetXmin()) / fPeriod );
		FinishSample    = floor( (((*Pulses)[i]).fTime - fDelay + fPMT->GetXmax()) / fPeriod );
		// Get delay time from "0" of OutWave to "0" of SPE shape
		PulseTime       = ((*Pulses)[i]).fTime - fDelay;
		// Get amplitude of SPE shape
		PulseAmpl       = ((*Pulses)[i]).fAmpl;
		// Limit edges
		if (StartSample < 0)
			StartSample = 0;
		if (FinishSample > fNumSamples - 1)
			FinishSample = fNumSamples - 1;
		// Add SPE to OutWave
		for (int s = StartSample; s <= FinishSample; s++) {
			SampleTime = s*fPeriod;
			fOutWave[s] += PulseAmpl/fGain * fPMT->Eval(SampleTime - PulseTime);
		}
	}
}

// Draw histograms
void MakeWave::DrawHists() {

	// CREATE HISTOGRAMS

	// For fPhotoElectrons
	// if (!TimeHist) {
	// 	Double_t LowTime  = 0; // Lower bound of time delay for histogram;
	// 	Double_t HighTime = fPMT->GetTOFe() + 3 * fPMT->GetTOFe_Sigma(); // Upper bound of time delay for histogram;
	// 	TimeHist = new TH1F ("TimeHist", "Time delay from photon hit to pulse", 100, LowTime, HighTime);
	// }
	// if (!AmplHist) {
	// 	Double_t LowAmpl  = 0; //GetAmpl() - 3 * GetAmpl_Sigma(); // Lower bound of amplitude for histogram;
	// 	Double_t HighAmpl = fPMT->GetAmpl() + 3 * fPMT->GetAmpl_Sigma(); // Upper bound of amplitude for histogram;
	// 	AmplHist = new TH1F ("AmplHist", "Amplitude of pulse", 100, LowAmpl, HighAmpl);
	// }
	// TCanvas *c4 = new TCanvas();
	// c4->SetTitle("Distribution of amplitude and time for Photon Pulses");
	// c4->Divide(2,1);
	// c4->cd(1);
	// R11->AmplHist->Draw();
	// c4->cd(2);
	// R11->TimeHist->Draw();

	// For pulse area
	if (!fPulseAreaHist) {
		Double_t Low_Area  = 0;
		Double_t High_Area = 3*fPMT->GetArea()/(ns*mV);
		fPulseAreaHist  = new TH1F ("fPulseAreaHist","PulseArea[mV*ns]",1000, Low_Area, High_Area);
	}
	for (unsigned int i=0; i< fPhotoElectrons->size(); i++) {
		RED::PMT::Pulse OnePulse = (*fPhotoElectrons).at(i);
		fPulseAreaHist->Fill (OnePulse.fAmpl * fPMT->GetShapeArea() / (mV*ns));
	}
	TCanvas *c3 = new TCanvas();
	c3->SetTitle("Pulse Area");
	c3->cd();
	c3->SetLogy();
	fPulseAreaHist->Draw();

	// For dark electrons
	if (!fDarkTimeHist) {
		Double_t LowDarkTime  = fDelay;
		Double_t HighDarkTime = fDelay + fPeriod*fNumSamples;
		fDarkTimeHist = new TH1F ("fDarkTimeHist", "Abs time of dark pulse", 100, LowDarkTime, HighDarkTime);
	}
	if (!fDarkAmplHist) {
		Double_t LowDarkAmpl  = fPMT->GetAmpl() - 3 * fPMT->GetAmpl_Sigma();
		Double_t HighDarkAmpl = fPMT->GetAmpl() + 3 * fPMT->GetAmpl_Sigma();
		fDarkAmplHist = new TH1F ("fDarkAmplHist", "Amplitude of dark pulse", 100, LowDarkAmpl, HighDarkAmpl);
	}
	for (unsigned int i =0; i< fDarkElectrons->size(); i++) {
		RED::PMT::Pulse DarkPulse = (*fDarkElectrons).at(i);
		fDarkTimeHist->Fill (DarkPulse.fTime);
		fDarkAmplHist->Fill (DarkPulse.fAmpl);
	}
	TCanvas *c5 = new TCanvas();
	c5->SetTitle("Distribution of amplitude and time for Dark Pulses");
	c5->Divide(2,1);
	c5->cd(1);
	fDarkAmplHist->Draw();
	c5->cd(2);
	fDarkTimeHist->Draw();
	c5->WaitPrimitive();
	
	// For number of fPhotoElectrons
	TCanvas *c2 = new TCanvas();
	c2->SetTitle("Number of phe");
	c2->cd();
	c2->SetLogy();
	fNumPheHist->Draw();
	fNumPheHist->SetTitle("Number of phe resulted (-1 or -2 == 1 or 2 phe from 1st dynode);N(phe);Events");
}

// Print OutWave
void MakeWave::PrintOutWave() {
	cout << "Printing OutWave..." << endl;
	for (int i = 0; i < fNumSamples; i++) {
		cout << fOutWave[i] << "\t(t = " << (fDelay + i*fPeriod)/ns << " ns)" << endl;
	}
	cout << "OutWave was printed" << endl;
}

// Draw OutWave
void MakeWave::DrawOutWave () {
	cout << "Drawing OutWave..." << endl;
	TGraph  *g1 = new TGraph(fNumSamples);
	for (int i = 0; i < fNumSamples; i++) {
		g1->SetPoint(i, (fDelay + i*fPeriod)/ns, fOutWave[i]);
	}
	g1->SetTitle("Waveform");
	g1->Draw();
	cout << "OutWave was drawn" << endl;
}

void MakeWave::SaveOutWave (const char *filename) {
	RED::OutputFile *outfile = new RED::OutputFile(filename);
	outfile->Open();
	if (outfile->IsOpen()) {
		RED::Event *event = new RED::Event;
		RED::Waveform *wf = event->GetNewWaveform();
		event->fNumChannels = event->GetNumWaveforms();
		event->SetRunInfo(new RED::RunInfo);
		wf->fChannel = 0;
		wf->fPeriod  = fPeriod;
		wf->fGain    = fGain;
		wf->fNumSamples = fNumSamples;
		wf->fDelay   = fDelay;
		wf->fData.assign(fOutWave.begin(), fOutWave.end());
		wf->Print();
		outfile->WriteEvent(event);
		cout << "ready" << endl;
	}
	else  {
		cout << "ERROR. File can't be writed" << endl;
	}
	outfile->Close();
}