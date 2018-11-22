#include <iostream>

#include <TRandom3.h>
#include "SystemOfUnits.h"
#include <TCanvas.h>
#include <TGraph.h>

#include "MakeTest.h"

using std::vector;
using std::cout;
using std::endl;
using CLHEP::ns;

MakeTest::MakeTest() {
	fTimesVec = new vector <Double_t>;
	fMeanOutWave = new vector <Double_t>;
	cout << "MakeTest object created" << endl;
}

void MakeTest::SetPhotonsTimes (Double_t *TimesArr) {
	Int_t TimesArrSize = sizeof(TimesArr)/sizeof(TimesArr[0]);
	(*fTimesVec).resize (TimesArrSize);
	for (int i = 0; i < TimesArrSize; i++) {
		(*fTimesVec)[i] = TimesArr[i]*ns;
	}
}

void MakeTest::AverageOW (Int_t DebugN, MakeWave* MakeWaveObj) {

	//Double_t Period     = MakeWaveObj -> GetPeriod();
	//Double_t Gain       = MakeWaveObj -> GetGain();
	Double_t NumSamples = MakeWaveObj -> GetNumSamples();
	//Double_t Delay      = MakeWaveObj -> GetDelay();

	cout << "Create vector for mean waveform" << endl;
	fMeanOutWave->resize (NumSamples);
	Bool_t ShowAllNum = false; // Report status of each OutWave performing
	                           // If false - only 10 steps will be showed
	// Create OutWave
	if (DebugN < 100) {
		ShowAllNum = true;
	}
	cout << "Get into cycle" << endl;
	for (int i = 1; i <= DebugN; i++){
		if (ShowAllNum)
			cout << "Creating " << i << " waveform" << endl;
		else{
			if (!(i % int(floor(DebugN/10)))){
				cout << "Creating " << i << " waveform" << endl;
			}
		}
		MakeWaveObj->CreateOutWave();
		for (int k = 0; k < NumSamples; k++){
			(*fMeanOutWave)[k] += MakeWaveObj->GetOutWave()[k];
		}
	}
}

void MakeTest::DrawOW (MakeWave* MakeWaveObj) {
	Double_t Period     = MakeWaveObj -> GetPeriod();
	//Double_t Gain       = MakeWaveObj -> GetGain();
	Double_t NumSamples = MakeWaveObj -> GetNumSamples();
	Double_t Delay      = MakeWaveObj -> GetDelay();

	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->SetTitle("OutWave");
	TGraph *g2  = new TGraph (NumSamples);
	for (int i = 0; i < NumSamples; i++){
		g2->SetPoint (i, (Delay + i*Period)/ns, (*fMeanOutWave)[i]);
		// cout << "point num " << i << "at time= " << (Delay + i*Period)/ns << " ns ";
		// cout << "with ampl= " << (*fMeanOutWave)[i] << " ADC-units was added" <<endl;
	}
	g2->SetTitle("OutWave;Time, ns;Amplitude, ADC units");
	g2->Draw();
}

void MakeTest::RandPhotonTimes (Int_t number, Double_t leftEdge, Double_t rightEdge) {
	fTimesVec -> resize (number);
	TRandom3 rand;
	rand.SetSeed(0);
	for (int i=0; i<number; i++){
		(*fTimesVec)[i] = leftEdge + rand.Rndm() * (rightEdge - leftEdge);
	}
	sort((*fTimesVec).begin(), (*fTimesVec).end());
}

void MakeTest::PrintTimesVec () {
	cout << "Arrival times:" << endl;
	for (unsigned int i = 0; i < fTimesVec->size(); i++){
		cout << ((*fTimesVec)[i])/ns << " ns" << endl;
	}
}

void MakeTest::DrawShape (RED::PMT *pmt, Char_t const *title) {
	TCanvas *c5 = new TCanvas();
	c5->SetTitle("SPE Shape");
	c5->cd();
	switch (pmt->fMode) {
		case RED::PMT::kModeNone :
			break;
		case RED::PMT::kModeF1 :
			pmt->fShape.func->SetTitle (title);
			pmt->fShape.func->Draw();
			break;
		case RED::PMT::kModeSpline :
			pmt->fShape.spline->SetTitle (title);
			pmt->fShape.spline->Draw();
			break;
	}
	//c5->WaitPrimitive();
}