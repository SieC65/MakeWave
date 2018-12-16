#include <iostream>

#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include "SystemOfUnits.h"

#include "MakeWave.h"

using std::cout;
using std::endl;
using CLHEP::ns;

// Simple constructor
MakeWave::MakeWave () {
	fPhotoElectrons   = 0;
	fDarkElectrons    = 0;
	fPMT              = 0;
	fPulseAreaHist    = 0;
	fOutFile          = 0;
}

// Set PMT
void MakeWave::SetPMT (RED::PMT* pmt) {
	fPMT = pmt;
}

// Set OutWave parameters
void MakeWave::SetOutWave (Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay) {
	fPeriod     = Period;
	fGain       = Gain;
	fNumSamples = NumSamples;
	fDelay      = Delay;
}

// Set some default parameters
void MakeWave::SetDefaults(){
	fPeriod     = 2*ns;
	fGain       = 0.125*mV;
	fNumSamples = 150000;
	fDelay      = -150000*ns;
}

// Set sequence of photon times
void MakeWave::SetPhotonTimes (vector <double>* PhotonTimes) {
	fPhotonTimes = PhotonTimes;
}

Double_t MakeWave::GetFrac (Double_t FracWindow, Double_t TotalWindow) {

	if (fPhotoElectrons->size()) {
		// Find min of delays
		Double_t minDelay = (fPhotoElectrons->at(0)).fTime;
		for (unsigned int i = 0; i < fPhotoElectrons->size(); i++) {
			if (minDelay > (fPhotoElectrons->at(i)).fTime)
				minDelay = (fPhotoElectrons->at(i)).fTime;
		}
		// Calculate integrals in both windows
		Double_t TotalSum = 0;
		Double_t FracSum = 0;
		if (!TotalWindow) {		
			for (unsigned int i = 0; i < fPhotoElectrons->size(); i++) {
				TotalSum += (fPhotoElectrons->at(i)).fAmpl;
				if ((fPhotoElectrons->at(i)).fTime - minDelay < FracWindow)
					FracSum += (fPhotoElectrons->at(i)).fAmpl;
			}
		}
		else {
			for (unsigned int i = 0; i < fPhotoElectrons->size(); i++) {
				if ((fPhotoElectrons->at(i)).fTime - minDelay < TotalWindow) {
					TotalSum += (fPhotoElectrons->at(i)).fAmpl;
					if ((fPhotoElectrons->at(i)).fTime - minDelay < FracWindow)
						FracSum += (fPhotoElectrons->at(i)).fAmpl;
				}
			}
		}
		Double_t frac = FracSum / TotalSum;
		return frac;
	}
	else
		return 0;
}

// Creating OutWave
void MakeWave::CreateOutWave () {
	fOutWave.clear ();                //  Clear vector OutWave
	fOutWave.resize (fNumSamples, 0); // Resize vector OutWave

	// ADD SPE FROM PHOTONS

	// Check if PulseArray vector fPhotoElectrons exists
	if (!fPhotoElectrons)
		fPhotoElectrons = new RED::PMT::PulseArray;
	else fPhotoElectrons->clear();

	// Generate fPhotoElectrons, add them to OutWave
	for (unsigned int i = 0; i < fPhotonTimes->size(); i++) {
		fPMT->OnePhoton (fPhotonTimes->at(i), *fPhotoElectrons, false);
	}
	AddPulseArray (fPhotoElectrons);

	// ADD DARK COUNTS

	// Check if PulseArray vector fDarkElectrons exists
	if (!fDarkElectrons)
		fDarkElectrons = new RED::PMT::PulseArray;
	else fDarkElectrons->clear();

	// Generate dark electrons and add them to OutWave
	fPMT->GenDCR (fDelay - (fPMT->GetXmax() - fPMT->GetXmin()), fDelay + fNumSamples * fPeriod, *fDarkElectrons);
	AddPulseArray (fDarkElectrons);

	//cout << "OutWave was created" << endl;
}

// Print OutWave
void MakeWave::PrintOutWave() {
	cout << "Printing OutWave..." << endl;
	for (int i = 0; i < fNumSamples; i++)
		cout << fOutWave[i] << "\t(t = " << (fDelay + i*fPeriod)/ns << " ns)" << endl;
	cout << "End of OutWave" << endl;
}

// Draw histograms
void MakeWave::DrawHists() {

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
}

// Draw OutWave
void MakeWave::DrawOutWave () {
	TCanvas *c1 = new TCanvas();
	c1->cd();
	TGraph  *g1 = new TGraph(fNumSamples);
	for (int i = 0; i < fNumSamples; i++) {
		g1->SetPoint(i, (fDelay + i*fPeriod)/ns, fOutWave[i]);
	}
	g1->SetTitle("Waveform");
	g1->Draw();
}

void MakeWave::AddToFile () {
	if (fOutFile) {
		if (fOutFile->IsOpen()) {
			fWaveform->fData.assign(fOutWave.begin(), fOutWave.end());
			fOutFile->WriteEvent(fEvent);
			fNumEv++;
		}
		else  {
			cout << "ERROR. File can't be written" << endl;
		}
	}
	else {
		cout << "ERROR. File was not created" << endl;
	}
}

RED::OutputFile* MakeWave::GetNewFile(const char *filename) {
	if (fOutFile) {
		if (fOutFile->IsOpen()) {
			fOutFile->Close();
		}
	}
	fOutFile = new RED::OutputFile(filename);
	fOutFile->Open();
	if (fOutFile->IsOpen()){
		fEvent = new RED::Event;
		fRunInfo = new RED::RunInfo;
		fEvent->SetRunInfo (fRunInfo);
		fWaveform = fEvent->GetNewWaveform();
		fWaveform->fChannel = 0;
		fWaveform->fPeriod  = fPeriod;
		fWaveform->fGain    = fGain;
		fWaveform->fNumSamples = fNumSamples;
		fWaveform->fDelay   = fDelay;
		fEvent->fNumChannels = fEvent->GetNumWaveforms();;
		fNumEv = 0;
		return fOutFile;
	}
	else  {
		cout << "ERROR. File " << filename << " can't be written" << endl;
		return 0;
	}
}

void MakeWave::CloseFile() {
	if (fOutFile) {
		if (fOutFile->IsOpen()) {
			RED::RunInfo *info = new RED::RunInfo();
			info->fDAQ_ID="CENNS10";
			fOutFile->AddRunInfo(info);
			fOutFile->Close();
		}
	}
}

// Add PulseArray vector to OutWave
void MakeWave::AddPulseArray (RED::PMT::PulseArray *Pulses) {
	Double_t SampleTime   = 0; // Time of sample from "0" of OutWave
	Double_t PulseTime    = 0; // Time from "0" of OutWave to "0" of SPE shape
	Double_t PulseAmpl    = 0; // Amplitude of SPE shape
	Int_t StartSample     = 0; // First sample in SPE domain
	Int_t FinishSample    = 0; // Last sample in SPE domain
	
	// Go along all pulses and add them to OutWave
	for (unsigned int i = 0; i < Pulses->size(); i++) {
		// Get delay time from "0" of OutWave to "0" of SPE shape
		PulseTime    = (Pulses->at(i)).fTime - fDelay;
		// Calculate left and right samples including SPE
		StartSample  =  ceil( (PulseTime + fPMT->GetXmin()) / fPeriod );
		FinishSample = floor( (PulseTime + fPMT->GetXmax()) / fPeriod );
		// Get amplitude of SPE shape
		PulseAmpl    = (Pulses->at(i)).fAmpl;
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