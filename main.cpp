#include <iostream>
#include "MakeWave.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

using namespace CLHEP;

//Print 2 graphs (new and old algorithms) for sum of DebugN OutWaves
void Compare (int DebugN, const MakeWave::OutWavePar &OW, MakeWave *aex);
//Struct for random set of timeseq - vector of photon arrival times
struct ArrTimes {
	int number;			//number of arr. times
	Double_t lrange;	//left range of times
	Double_t rrange;	//right range of times
};
//Set random timeseq
void SetRandTimes (vector <double> *Tseq, ArrTimes AT);

int main() {
	Int_t DebugAlg		= 1000;	//0 - for only new algorithm, n - for Compare (n)
	Bool_t RandomTimes	= 1;	//0 - for manual, 1 - for random set photon times
	ArrTimes AT;	//Struct for random set of timeseq
		AT.number = 100;		//number of arr. times
		AT.lrange = -1000*ns;	//left range of times
		AT.rrange = 1000*ns;		//right range of times
	//Set SPE and OutWave parameters
	MakeWave::SPEPar SPEX;		//SPE Example
		//Parameters of SPE form
		SPEX.Type 			= 0;			//0 - gaussian , 1 - square pulse
		SPEX.Trig	 		= 0.01;			//Coeff: Trigger for SPE = Trig*Ampl
		SPEX.Width 			= 10*ns;		//FWHM of SPE gaussian (or width for square pulse)
		SPEX.Ampl 			= -1*mV;		//Amplitude of SPE
		SPEX.AmplSigma		= 0.2*mV;		//Sigma for SPE amplitude
		SPEX.Delay			= 20*ns;		//Delay of SPE signal
		SPEX.DelaySigma		= 5*ns;			//Sigma for SPE delay
		SPEX.Domain			= SPEX.Width;	//Domain width of SPE (initializing)
		//Case of real photon interaction with PMT
		SPEX.QE				= 0.3;		//Quantum Efficiency (full)
		SPEX.DPE_PC			= 0.225;	//Double Photoelectron Emission = P(2phe)/(P(2phe)+P(1phe)) for PC
		SPEX.QE_1d			= 0.105;	//Quantum Efficiency for 1dyn (only 1phe)
		SPEX.GF_1d			= 0.1;		//Arb. geometrical factor of 1dyn
		SPEX.ArbCurr_1d		= 0.2;		//Current(1d)/Current(PC) for SPE
		SPEX.ArbAmpl_1d		= SPEX.ArbCurr_1d;	//Ampl(1d)/Ampl(PC) for SPE
		SPEX.AmplSigma_1d	= SPEX.AmplSigma*SPEX.ArbAmpl_1d;	//Sigma for amplitude from 1dyn
		SPEX.Pphe_1d		= SPEX.GF_1d * SPEX.QE_1d;	//Probability of photoeffect on 1dyn for photons passed PC
		SPEX.QE_1d_ratio	= SPEX.ArbCurr_1d * SPEX.GF_1d * SPEX.QE_1d;	//auxiliary quantity
		SPEX.Pphe_PC		= (SPEX.QE - SPEX.QE_1d_ratio)/(1 + SPEX.DPE_PC - SPEX.QE_1d_ratio);	//Probability of photoeffect on PC
		SPEX.Ampl_DPE_PC	= 2 * SPEX.Ampl;	//Amplitude of signal from 2phe from PC
	MakeWave::OutWavePar OWX;	//OutWave example
		OWX.Period 		= 0.2*ns;	//Time between samples of OutWave
		OWX.Gain 		= 0.125*mV;	//Units of ADC
		OWX.Num 		= (AT.rrange-AT.lrange)/OWX.Period;		//Number of samples in OutWave
		OWX.Delay 		= AT.lrange*ns;		//Delay from "0" of abs.time to "0" sample of OutWave

	vector <double> Tseq;
	if (!RandomTimes) {
		//Set sequence of SPE arrival times
		Double_t TseqArr[] = { -20 , 0 , 5.5 , 15 , 50};
		Int_t TseqArrSize = sizeof(TseqArr) / sizeof(TseqArr[0]);
		Tseq.resize (TseqArrSize);
		for (int i = 0; i < TseqArrSize; i++) {
			TseqArr[i]	*= ns;	//Set times in ns
			Tseq[i]		= TseqArr[i];
		}
	}
	else {
		Tseq.resize (AT.number);
		SetRandTimes(&Tseq, AT);
	}
//	cout << "Times of photons:" << endl;
	for (int i=0; i< Tseq.size(); i++) {
	//	cout << Tseq[i]/ns << " ns" << endl;
	}

	//Set parameters given above and create OutWave
	MakeWave *a = new MakeWave;	//Create object of MakeWave class
	TApplication *app = new TApplication("Canvas",0,0);
	a->SetSPE (SPEX);			//Set SPE
	a->SetParams (OWX);		//Set OutWave parameters
	a->SetTimeSeq (&Tseq);		//Set sequence of photon arrival times
	a->SetRand(true);			//Randomize random generator for SPE params
	if (DebugAlg == 0) {
		a->CreateOutWave();		//Just create OutWave vector
		TCanvas *c3 = new TCanvas();
		c3->cd();
		a->NPhe->Draw();
		TCanvas *c4 = new TCanvas();
		c4->cd();
		c4->SetLogy();
		a->PA->Draw();
		a->Draw();				//And draw OutWave on screen
	}
	else {
		Compare (DebugAlg, OWX, a);	//Debug algorithm
	}
	cout << "It was " << (DebugAlg + !DebugAlg)*Tseq.size();
	cout << " photons that created:" << endl;
	cout << "0 phe in PC:" << a->NPhe->GetBinContent(1) << endl;
	cout << "1 phe in PC:" << a->NPhe->GetBinContent(2) << endl;
	cout << "2 phe in PC:" << a->NPhe->GetBinContent(3) << endl;
	cout << "1 phe in 1d:" << a->NPhe->GetBinContent(4) << endl;
	return 0;
}

//DEBUG FUNC
//Create average OutWave for DebugN OutWaves with the same SPE arrival times
void Compare (int DebugN, const MakeWave::OutWavePar &OW, MakeWave *aex) {
	vector <double> MeanOutWaveNew;		//sum of DebugN OutWaves from new alg.
	vector <double> MeanOutWaveOld;		//sum of DebugN OutWaves from old alg.
	MeanOutWaveNew.resize(OW.Num);
	MeanOutWaveOld.resize(OW.Num);
	for (int i = 1; i <= DebugN; i++) {
		cout << "Creating " << i << " waveform" << endl;
		//For old algorithm
		aex->CreateOutWaveOld();
		for (int k = 0; k < OW.Num; k++) {
			MeanOutWaveOld[k] += aex->OutWave[k];
		}
		//For new algorithm
		aex->CreateOutWave();
		for (int k = 0; k < OW.Num; k++) {
			MeanOutWaveNew[k] += aex->OutWave[k];
		}
	}
	TCanvas *c2 = new TCanvas();
	c2->Divide(1,2);
	//Draw OutWave for new algorithm
	c2->cd(1);
	TGraph *g2 = new TGraph(OW.Num);
	for (int i = 0; i < OW.Num; i++)
		g2->SetPoint(i, (OW.Delay + i*OW.Period)/ns, MeanOutWaveNew[i]);
	g2->SetTitle("New");
	g2->Draw();
	//Draw OutWave for old algorithm
	c2->cd(2);
	TGraph *g3 = new TGraph(OW.Num);
	for (int i = 0; i < OW.Num; i++)
		g3->SetPoint(i, (OW.Delay + i*OW.Period)/ns, MeanOutWaveOld[i]);
	g3->SetTitle("Old");
	g3->Draw();	
	TCanvas *c3 = new TCanvas();
	c3->cd();
	aex->NPhe->Draw();
	TCanvas *c4 = new TCanvas();
	c4->cd();
	c4->SetLogy();
	aex->PA->Draw();
	c4->WaitPrimitive();
}

void SetRandTimes (vector <double> *Tseq, ArrTimes AT) {
	TRandom3 rand;
	rand.SetSeed(0);
	for (int i=0; i<AT.number; i++) {
		(*Tseq)[i] = AT.lrange + rand.Rndm()*(AT.rrange-AT.lrange);
	}
}