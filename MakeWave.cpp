#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include "SystemOfUnits.h"

using namespace CLHEP;

MakeWave::MakeWave () : fisrnd(true) {
	cout << "MakeWave object was created" << endl;
}

void MakeWave::SetSPE (struct SPEPar SPEX) {
	fSPEX = SPEX;
	//SPE form is set below
	if (fSPEX.Type == 1) {
		//Square pulse
		fSPEX.Domain = fSPEX.Width;				//Domain of SPE
		fSPE = new TF1("SPE","((x > [0]) && (x < [1]))*[2]");
		fSPE->SetParameter(0, 0);				//Left end of the range
		fSPE->SetParameter(1, fSPEX.Domain);	//Right end of the range
		fSPE->SetParameter(2, fSPEX.Ampl);		//Amplitude of square pulse
	}
	else {
		//Gaussian
		//Domain of SPE (region with value > fSPEX.Trig*fSPEX.Ampl)
		fSPEX.Domain = fSPEX.Width * (sqrt(log2(1/fSPEX.Trig)) + 1)/2;	
		fSPE = new TF1("SPE","gaus(0)",0,fSPEX.Domain);
		fSPE->SetParameter(0, fSPEX.Ampl);		//Amplitude of gauss
		fSPE->SetParameter(1, fSPEX.Width/2);	//Center
		fSPE->SetParameter(2, (fSPEX.Width)/(2*sqrt(2*log(2))));	//Sigma
		fSPE->Draw();
	}	
	fSPE->SetNpx(1000);
	cout << "SPE signal was set" << endl;
}

void MakeWave::SetParams (struct OutWavePar OWX) {
	fOWX = OWX;
	cout << "Parameters were set" << endl;
}

void MakeWave::SetTimeSeq (vector <double> *Tseq) {
	ftimeseq = Tseq;	//now addresses of ftimeseq and Tseq are equal (pointer Tseq created in main.cpp)
	cout << "Sequence of SPE's arrival times was set" << endl;
}

void MakeWave::SetRand(Bool_t isrnd) {
	fisrnd = isrnd;
	if (fisrnd == true)
		fRND.SetSeed(0);
}

Bool_t MakeWave::GetIsrnd () {
	return fisrnd;
}

void MakeWave::CreateOutWave() {
	Double_t sampletime = 0;
	OutWave.clear();
	OutWave.resize (fOWX.Num, 0);
	Double_t SPEAmplSample	= 0;
	Double_t SPEDelaySample	= 0;
	Int_t startsample		= 0;
	Int_t finishsample		= 0;
	Bool_t PosAmplSigma		= true;
	Bool_t PosDelaySigma	= true;
	if (fSPEX.AmplSigma <= 0) {
		PosAmplSigma	= false;
		SPEAmplSample	= 1;
	}
	if (fSPEX.DelaySigma <= 0) {
		PosDelaySigma	= false;
		SPEDelaySample	= fSPEX.Delay;
	}
	cout << "SPE Domain=" << fSPEX.Domain << endl;
	for (int i = 0; i < int(ftimeseq->size()); i++) {
		if (PosAmplSigma)
			SPEAmplSample	= fRND.Gaus(1, fSPEX.AmplSigma/fSPEX.Ampl);
		if (PosDelaySigma)
			SPEDelaySample	= fRND.Gaus(fSPEX.Delay, fSPEX.DelaySigma);
		startsample 	= int(((*ftimeseq)[i] + SPEDelaySample - fOWX.Delay)/fOWX.Period) + 1;	//first sample in SPE domain
		finishsample 	= int(((*ftimeseq)[i] + SPEDelaySample - fOWX.Delay + fSPEX.Domain)/fOWX.Period);	//last sample in SPE domain
		if (startsample < 0)
			startsample = 0;
		if (finishsample > fOWX.Num - 1)
			finishsample = fOWX.Num - 1;
		cout << "for " << i << " SPE AMPL=" << SPEAmplSample << ", DELAY=" << SPEDelaySample;
		cout << ", ARRTIME=" << (*ftimeseq)[i] + SPEDelaySample - fOWX.Delay << endl;		
		cout << "Write from " << startsample << " to " << finishsample << " samples";
		cout << " (Domain = " << fSPEX.Domain << ")" << endl;
		for (int sample = startsample; sample <= finishsample; sample++) {
			sampletime = fOWX.Delay + sample*fOWX.Period;
			OutWave[sample] += SPEAmplSample*(fSPE->Eval(sampletime - SPEDelaySample - (*ftimeseq)[i]))/fOWX.Gain;
			cout << "\tfor " << sample << " sample (" << sampletime/ns << " ns) added ";
			cout << SPEAmplSample*(fSPE->Eval(sampletime - SPEDelaySample - (*ftimeseq)[i]))/fOWX.Gain << endl;
		}
	}
//	cout << "OutWave was created" << endl;
}

void MakeWave::CreateOutWaveOld() {
	Double_t t = 0;
	OutWave.clear();
	OutWave.resize (fOWX.Num, 0);
	vector <double> SPEAmplVec(ftimeseq->size());
	vector <double> SPEDelayVec(ftimeseq->size());	
	for (int i = 0; i < int(ftimeseq->size()); i++) {
	//	cout << "before random" << endl;
	//	cout << "for " << i << " SPE AMPL=" << SPEAmplVec[i] << "and DEL=" << SPEDelayVec[i] << endl;
		if (fisrnd == true) {
			fRND.SetSeed(0);
		}
		SPEAmplVec[i] 	= fRND.Gaus(1, fSPEX.AmplSigma/fSPEX.Ampl);
		SPEDelayVec[i] 	= fRND.Gaus(fSPEX.Delay, fSPEX.DelaySigma);
	//	cout << "after random" << endl;
	//	cout << "for " << i << " SPE AMPL=" << SPEAmplVec[i] << "and DEL=" << SPEDelayVec[i] << endl;
	}
	for (int sample = 0; sample < fOWX.Num; sample++) {
		t = fOWX.Delay + sample*fOWX.Period;
		for (int i = 0; i < int(ftimeseq->size()); i++) {
			OutWave[sample] += (SPEAmplVec[i])*(fSPE->Eval(t - SPEDelayVec[i] - (*ftimeseq)[i]))/fOWX.Gain;
		//	cout << "for " << sample << " sample added ";
		//	cout << (SPEAmplVec[i])*(fSPE->Eval(t - SPEDelayVec[i] - (*ftimeseq)[i]))/fOWX.Gain << endl;
		}
	}
//	cout << "OutWave was created"<< endl;
}

void MakeWave::PrintOutWave() {
	for (int i = 0; i < fOWX.Num; i++) {
		cout << OutWave[i] << "\t(t = " << (fOWX.Delay + i*fOWX.Period)/ns << " ns)" << endl;
	}
	cout <<"OutWave was printed"<< endl;
}

void MakeWave::Draw (TString name) {
	TCanvas *c1 = new TCanvas();
	TGraph *g1 = new TGraph(fOWX.Num);
	for (int i = 0; i < fOWX.Num; i++) {
		g1->SetPoint(i, (fOWX.Delay + i*fOWX.Period)/ns, OutWave[i]);
	}
	g1->Draw();
	if (name != "")
		c1->SaveAs(name);	//graph is closed when enabled
	c1->WaitPrimitive();
}