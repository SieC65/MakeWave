#include <iostream>
#include "MakeWave.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

using namespace CLHEP;

//Print 2 graphs (new and old algorithms) for sum of DebugN OutWaves
void Compare (int DebugN, const MakeWave::OutWavePar &OW, MakeWave *aex);
//Struct for setting timeseq - vector of photon arrival times
struct ArrTimes {
	int number;			//Number of arr. times. Manual setting if =0
	Double_t lrange;	//Left range of times
	Double_t rrange;	//Right range of times
};
//Set timeseq
void SetTseq (ArrTimes AT, vector <double> *Tseq, Double_t *TseqArr, Int_t TseqArrSize);

int main() {
	Int_t DebugAlg		= 0;	//0 - for only new algorithm, n - for Compare (n)
	ArrTimes AT;	//Struct for setting timeseq
		AT.number = 100;	//Number of arrival times: 0 - for manual set times, N - for random set
		AT.lrange = -1000*ns;	//Left range of times
		AT.rrange =  1000*ns;	//Right range of times
	//Set SPE and OutWave parameters
	MakeWave::SPEPar SPEX;			//SPE Parameters
		//Parameters of SPE form
		SPEX.Type 		= 0;			//0 - gaussian , 1 - square pulse
		SPEX.Trig	 	= 0.01;			//"Trigger": ampl of SPE at arr. times = Trig*Ampl
		SPEX.Width 		= 10*ns;		//FWHM of SPE gaussian (or width for square pulse)
		SPEX.Ampl 		= -1*mV;		//Amplitude of SPE
		SPEX.AmplSigma	= 0.2*mV;		//Sigma for SPE amplitude
		SPEX.Delay		= 20*ns;		//Delay of SPE signal
		SPEX.DelaySigma	= 0*ns;			//Sigma for SPE delay
		SPEX.Domain		= SPEX.Width;	//Domain width of SPE (initializing)
	MakeWave::OutWavePar OWX;		//OutWave Parameters
		OWX.Period 		= 2*ns;			//Time between samples of OutWave
		OWX.Gain 		= 0.125*mV;		//Units of ADC
		OWX.Num 		= ceil((AT.rrange-AT.lrange)/OWX.Period);	//Number of samples in OutWave
		OWX.Delay 		= AT.lrange*ns + SPEX.Delay;	//Delay from "0" of abs.time to "0" sample of OutWave
	
	//Set PMT Hamamatsu R11410-20 parameters
	//Some values are taken from articles
	//C.H.Faham et.al. "Measurements of wavelength-dependent double photoelectron emission..." //JINST 2015
	//and
	//D.Yu.Akimov et.al. "Performance of Hamamatsu R11410-20 PMTs..." //JINST 2016
	PMT *R11 = new PMT;
	PMT::PMTPar PMTX;	//PMT Parameters
		PMTX.QE				= 0.3;		//Quantum Efficiency (full)
		PMTX.DPE			= 0.225;	//Double Photoelectron Emission = P(2phe)/(P(2phe)+P(1phe)) for PC
		PMTX.QE_1d			= 0.105;	//Quantum Efficiency for 1dyn (only 1phe)
		PMTX.ArbAmpl_1d		= (1./13);		//Ampl(1d)/Ampl(PC) for SPE
		PMTX.GF_1d			= 0.1;		//Arb. geometrical factor of 1dyn
		PMTX.TOF_1d			= 0.2*ns;	//TOF for photon from PC to 1dyn
	R11->SetPMTPar (PMTX, SPEX.AmplSigma);

	//Create Timeseq
	vector <double> Tseq;
	Double_t TseqArr[] = { -20 , 0 , 5.5 , 15 , 50};	//Manual values
	SetTseq (AT, &Tseq, TseqArr, sizeof(TseqArr)/sizeof(Double_t));
	
	//Set parameters given above and create OutWave
	MakeWave *a = new MakeWave;	//Create object of MakeWave class
	TApplication *app = new TApplication("Canvas",0,0);
	a->SetSPE (SPEX);			//Set SPE
	a->SetPMT (R11);			//Set used PMT
	a->SetOutWave (OWX);		//Set OutWave parameters
	a->SetTimeSeq (&Tseq);		//Set sequence of photon arrival times
	a->SetRand(true);			//Randomize random generator for SPE params

	if (DebugAlg == 0) {
	//Use only new algorithm
		a->CreateOutWave();		//Just create OutWave vector
		a->Draw();				//And draw OutWave on screen
	}
	else {
	//Debug algorithm
		Compare (DebugAlg, OWX, a);
	}
	cout << "It were " << (DebugAlg + !DebugAlg)*Tseq.size() << " photons that created:" << endl;
	cout << "0 phe in PC:" << a->NPhe->GetBinContent(1) << endl;
	cout << "1 phe in PC:" << a->NPhe->GetBinContent(2) << endl;
	cout << "2 phe in PC:" << a->NPhe->GetBinContent(3) << endl;
	cout << "1 phe in 1d:" << a->NPhe->GetBinContent(4) << endl;
	TCanvas *c3 = new TCanvas();
	c3->cd();
	a->NPhe->Draw();
	TCanvas *c4 = new TCanvas();
	c4->cd();
	c4->SetLogy();
	a->PA->Draw();
	c4->WaitPrimitive();
	return 0;
}

//DEBUG FUNC
//Create average OutWave for DebugN OutWaves with the same SPE arrival times
void Compare (int DebugN, const MakeWave::OutWavePar &OW, MakeWave *aex) {
	vector <double> MeanOutWaveNew;		//Sum of DebugN OutWaves from new alg.
	vector <double> MeanOutWaveOld;		//Sum of DebugN OutWaves from old alg.
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
}

void SetTseq (ArrTimes AT, vector <double> *Tseq, Double_t *TseqArr, Int_t TseqArrSize) {
	if (!AT.number) {
		//Manual set sequence of SPE arrival times
		cout << TseqArrSize << endl;
		(*Tseq).resize (TseqArrSize);
		for (int i = 0; i < TseqArrSize; i++) {
			(*Tseq)[i] = TseqArr[i]*ns;
		}
	}
	else {
		(*Tseq).resize (AT.number);
		TRandom3 rand;
		rand.SetSeed(0);
		for (int i=0; i<AT.number; i++) {
			(*Tseq)[i] = AT.lrange + rand.Rndm()*(AT.rrange-AT.lrange);
		}
	}
}