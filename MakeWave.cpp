#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <time.h>

MakeWave::MakeWave () : fisrnd(true) {
	cout << "MakeWave object was created" << endl;
}

void MakeWave::SetSPE (TF1 *SPE, Double_t SPEAmpl, Double_t SPEAmplSigma, Double_t SPEDelay, Double_t SPEDelaySigma, Double_t SPEDomain) {
	fSPE			= SPE;	//now addresses of fSPE and SPE are equal (pointer SPE created in main.cpp)	
	fSPEAmpl		= SPEAmpl;
	fSPEAmplSigma	= SPEAmplSigma;
	fSPEDelay		= SPEDelay;
	fSPEDelaySigma	= SPEDelaySigma;
	fSPEDomain		= SPEDomain;
	cout << "SPE signal was set" << endl;
}

void MakeWave::SetParams (Int_t Num, Double_t Delay, Double_t Period, Double_t Gain) {
	fPeriod		= Period;
	fNum		= Num;
	fDelay		= Delay;
	fGain		= Gain;
	cout << "Parameters were set" << endl;
}

void MakeWave::SetTimeSeq (vector <double> *Tseq) {
	ftimeseq = Tseq;	//now addresses of ftimeseq and Tseq are equal (pointer Tseq created in main.cpp)
	cout << "Sequence of SPE's arrival times was set" << endl;
}

void MakeWave::CreateOutWave() {
	Double_t sampletime = 0;
	OutWave.clear();
	OutWave.resize (fNum, 0);
	Double_t SPEAmplSample	= 0;
	Double_t SPEDelaySample	= 0;
	Int_t startsample		= 0;
	Int_t finishsample		= 0;
	if (fisrnd == true)
		fRND.SetSeed(0);
//	cout << "SPE Domain=" << fSPEDomain << endl;
	for (int i = 0; i < int(ftimeseq->size()); i++) {
		SPEAmplSample	= fRND.Gaus(1, fSPEAmplSigma/fSPEAmpl);
		SPEDelaySample	= fRND.Gaus(fSPEDelay, fSPEDelaySigma);
		startsample 	= int(((*ftimeseq)[i] + SPEDelaySample - fDelay)/fPeriod) + 1;	//first sample in SPE domain
		finishsample 	= startsample + int(fSPEDomain/fPeriod);
		if (startsample < 0)
			startsample = 0;
		if (finishsample > fNum - 1)
			finishsample = fNum - 1;
	//	cout << "for " << i << " SPE AMPL=" << SPEAmplSample << ", DELAY=" << SPEDelaySample;
	//	cout << ", ARRTIME=" << (*ftimeseq)[i] + SPEDelaySample - fDelay << endl;		
	//	cout << "Write from " << startsample << " to " << finishsample << " samples" << endl;
		for (int sample = startsample; sample <= finishsample; sample++) {
			sampletime = fDelay + sample*fPeriod;
			OutWave[sample] += SPEAmplSample*(fSPE->Eval(sampletime - SPEDelaySample - (*ftimeseq)[i]))/fGain;
		//	cout << "\tfor " << sample << " sample (" << sampletime/ns << " ns) added ";
		//	cout << SPEAmplSample*(fSPE->Eval(sampletime - SPEDelaySample - (*ftimeseq)[i]))/fGain << endl;
		}
	}
//	cout << "OutWave was created" << endl;
}

void MakeWave::CreateOutWaveOld() {
	Double_t t = 0;
	OutWave.clear();
	OutWave.resize (fNum, 0);
	vector <double> SPEAmplVec(ftimeseq->size());
	vector <double> SPEDelayVec(ftimeseq->size());	
	for (int i = 0; i < int(ftimeseq->size()); i++) {
	//	cout << "before random" << endl;
	//	cout << "for " << i << " SPE AMPL=" << SPEAmplVec[i] << "and DEL=" << SPEDelayVec[i] << endl;
		if (fisrnd == true) {
			fRND.SetSeed(0);
		}
		SPEAmplVec[i] 	= fRND.Gaus(1, fSPEAmplSigma/fSPEAmpl);
		SPEDelayVec[i] 	= fRND.Gaus(fSPEDelay, fSPEDelaySigma);
	//	cout << "after random" << endl;
	//	cout << "for " << i << " SPE AMPL=" << SPEAmplVec[i] << "and DEL=" << SPEDelayVec[i] << endl;
	}
	for (int sample = 0; sample < fNum; sample++) {
		t = fDelay + sample*fPeriod;
		for (int i = 0; i < int(ftimeseq->size()); i++) {
			OutWave[sample] += (SPEAmplVec[i])*(fSPE->Eval(t - SPEDelayVec[i] - (*ftimeseq)[i]))/fGain;
		//	cout << "for " << sample << " sample added ";
		//	cout << (SPEAmplVec[i])*(fSPE->Eval(t - SPEDelayVec[i] - (*ftimeseq)[i]))/fGain << endl;
		}
	}
//	cout << "OutWave was created"<< endl;
}

void MakeWave::PrintOutWave() {
	for (int i = 0; i < fNum; i++) {
		cout << OutWave[i] << "\t(t = " << (fDelay + i*fPeriod)/ns << " ns)" << endl;
	}
	cout <<"OutWave was printed"<< endl;
}

void MakeWave::Draw (TString name) {
	TCanvas *c1 = new TCanvas();
	TGraph *g1 = new TGraph(fNum);
	for (int i = 0; i < fNum; i++) {
		g1->SetPoint(i, (fDelay + i*fPeriod)/ns, OutWave[i]);
	}
	g1->Draw();
//	if (name != "")
//		c1->SaveAs(name);	//graph is closed if uncommented
	c1->WaitPrimitive();
}