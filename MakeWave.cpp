#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>

MakeWave::MakeWave () {
	cout << "MakeWave object was created" << endl;
}

void MakeWave::SetSPE (TF1 *SPE) {
	fSPE = SPE;	//now addresses of fSPE and SPE are equal (pointer SPE created in main.cpp)
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
	Double_t t = 0;
	OutWave.resize (fNum);
	for (int sample = 0; sample < fNum; sample++) {
		t = fDelay + sample*fPeriod;
		for (int i = 0; i < int(ftimeseq->size()); i++) {
			OutWave[sample] += (fSPE->Eval((t - ((*ftimeseq)[i]))/ns))/fGain;
			//cout << "for " << sample << " sample added " << fSPE->Eval(t - ((*ftimeseq)[i])) << endl;
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