#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>


MakeWave::MakeWave ():fSPEAmplDistr("SPEAmplDistr","gaus(0)",-1,3), fSPEDelayDistr("SPEDelayDistr","gaus(0)", 0, 1) {
	cout << "MakeWave object was created" << endl;
}

void MakeWave::SetSPE (TF1 *SPE, Double_t SPEAmpl, Double_t SPEAmplSigma, Double_t SPEDelay, Double_t SPEDelaySigma) {
	fSPE = SPE;	//now addresses of fSPE and SPE are equal (pointer SPE created in main.cpp)	
	cout << "SPE signal was set" << endl;
	
	fSPEAmplDistr.SetParameter(0,1);
	fSPEAmplDistr.SetParameter(1,1);
	fSPEAmplDistr.SetParameter(2,(SPEAmplSigma/SPEAmpl));
	fSPEAmplDistr.SetNpx(1000);
	fSPEDelayDistr.SetParameter(0,1);
	fSPEDelayDistr.SetParameter(1,SPEDelay);
	fSPEDelayDistr.SetParameter(2,SPEDelaySigma);
	fSPEDelayDistr.SetRange(SPEDelay - 3*SPEDelaySigma, SPEDelay + 3*SPEDelaySigma);
	fSPEDelayDistr.SetNpx(1000);
	cout << "Distributions for SPE signal was set" << endl;	
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
	Double_t t = 0;
	OutWave.clear();
	OutWave.resize (fNum, 0);
	vector <double> SPEAmplVec(ftimeseq->size());
	vector <double> SPEDelayVec(ftimeseq->size());
	
	for (int i = 0; i < int(ftimeseq->size()); i++) {
		cout << "before random" << endl;
		cout << "for " << i << " SPE AMPL=" << SPEAmplVec[i] << "and DEL=" << SPEDelayVec[i] << endl;
		TRandom3 rnd(0);
		rnd.SetSeed(ULong_t(time(0)));
		cout << rnd.Rndm() << endl;
		cout << rnd.GetSeed() << endl;
		cout << time(0) << endl;
		fSPEAmplDistr.GetRandom();
		fSPEDelayDistr.GetRandom();
		SPEAmplVec[i] 	= fSPEAmplDistr.GetRandom();
		SPEDelayVec[i] 	= fSPEDelayDistr.GetRandom();
		cout << "after random" << endl;
		cout << "for " << i << " SPE AMPL=" << SPEAmplVec[i] << "and DEL=" << SPEDelayVec[i] << endl;
	}
	for (int sample = 0; sample < fNum; sample++) {
		t = fDelay + sample*fPeriod;
		for (int i = 0; i < int(ftimeseq->size()); i++) {
			OutWave[sample] += (SPEAmplVec[i])*(fSPE->Eval(t - SPEDelayVec[i] - (*ftimeseq)[i]))/fGain;
		//	cout << "for " << sample << " sample added ";
		//	cout << (SPEAmplVec[i])*(fSPE->Eval(t - SPEDelayVec[i] - (*ftimeseq)[i]))/fGain << endl;
		}
	}
	cout << "OutWave was created"<< endl;
}

void MakeWave::PrintOutWave() {
	for (int i = 0; i < fNum; i++) {
		cout << OutWave[i] << "\t(t = " << (fDelay + i*fPeriod)/ns << " ns)" << endl;
	}
	cout <<"OutWave was printed"<< endl;
}

void MakeWave::Draw() {
	TCanvas *c1 = new TCanvas();
	TGraph *g1 = new TGraph(fNum);
	for (int i = 0; i < fNum; i++) {
		g1->SetPoint(i, (fDelay + i*fPeriod)/ns, OutWave[i]);
	}
	g1->Draw();
	c1->WaitPrimitive();
}