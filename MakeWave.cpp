#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>

using namespace CLHEP;

//Simple constructor
MakeWave::MakeWave () : fisrnd(true) {
	cout << "MakeWave object was created" << endl;
}

//Set SPE parameters and create histograms for PulseArea and Number of photons
void MakeWave::SetSPE (const SPEPar &SPEX) {
	fSPEX	= SPEX;		//Set SPE parameters
	
	//Set SPE form
	if (fSPEX.Type == 1) {
		//Square pulse
		fSPEX.Domain = fSPEX.Width;				//Domain of SPE
		fSPE = new TF1("SPE","((x > [0]) && (x < [1]))*[2]");
		fSPE->SetParameter(0, 0);				// Left end of SPE domain range
		fSPE->SetParameter(1, fSPEX.Domain);	//Right end of SPE domain range
		fSPE->SetParameter(2, fSPEX.Ampl);		//Amplitude of square pulse
	}
	else {
		//Gaussian (accessory, with amplitude = 1)
		//Domain of SPE (length of region with value > fSPEX.Trig*fSPEX.Ampl)
		fSPEX.Domain = fSPEX.Width * sqrt(log2(1/fSPEX.Trig));
		Double_t lrange = (fSPEX.Width - fSPEX.Domain)/2;	// Left end of SPE domain range
		Double_t rrange = (fSPEX.Width + fSPEX.Domain)/2;	//Right end of SPE domain range
	//	cout << "lrange = " << lrange << " and rrange = " << rrange << endl;
		fSPE = new TF1("SPE","gaus(0)",lrange,rrange);
		fSPE->SetParameter(0, 1);		//Amplitude of gaussian ( = 1 )
		fSPE->SetParameter(1, fSPEX.Width/2);	//Center
		fSPE->SetParameter(2, (fSPEX.Width)/(2*sqrt(2*log(2))));	//Sigma
		fSPE->Draw();
	}
	
	//Right end of PulseArea hist
	Double_t XmaxPA = 2 * fSPEX.Width * (abs(fSPEX.Ampl) + 3*fSPEX.AmplSigma)/(ns*mV);
	PA		= new TH1F("PA","PulseArea[mV*ns]",100,0,XmaxPA);	//Area under pulse SPE (DPE)
	NPhe	= new TH1F("NPhe","NumberOfPhe",4,0,4);			//Number of photons created 0, 1 and 2 phe
	cout << "SPE signal was set" << endl;
}

//Set OutWave parameters
void MakeWave::SetParams (const OutWavePar &OWX) {
	fOWX = OWX;
	cout << "Parameters were set" << endl;
}

//Set sequence of photon times
void MakeWave::SetTimeSeq (vector <double> *Tseq) {
	ftimeseq = Tseq;	//now addresses of ftimeseq and Tseq are equal (pointer Tseq created in main.cpp)
	cout << "Sequence of SPE's arrival times was set" << endl;
}

//Randomize simulating amplitudes and delays
void MakeWave::SetRand(Bool_t isrnd) {
	fisrnd = isrnd;
	if (fisrnd == true)
		fRND.SetSeed(0);
}

//Return status of randomizing
Bool_t MakeWave::GetIsrnd () {
	return fisrnd;
}

//New version for creating OutWave
void MakeWave::CreateOutWave () {
	
	//Initializing variables
	Double_t sampletime = 0;	//time of sample from "0" of OutWave
	Double_t arrtime	= 0;	//arrival time from "0" of OutWave
	OutWave.clear();				//clear vector OutWave
	OutWave.resize (fOWX.Num, 0);	//clear vector OutWave
	Double_t SPEAmplSample	= 0;	//random SPE amplitude
	Double_t SPEDelaySample	= 0;	//random SPE delay
	Int_t startsample		= 0;	//first sample in SPE domain
	Int_t finishsample		= 0;	// last sample in SPE domain
	Bool_t PosAmplSigma		= true;	//is sigma for amplitude >0
	Bool_t PosDelaySigma	= true;	//is sigma for delay >0
	Double_t Ampl			= 0;	//amplitude of signal from photon
	Double_t AmplSigma		= 0;	//sigma of amplitude
	Bool_t Interaction		= true;	//interaction between photon and PMT
	
	//If sigma for ampl or delay <= 0 , don't simulate
	if (fSPEX.AmplSigma <= 0) {
		PosAmplSigma	= false;
		SPEAmplSample	= fSPEX.Ampl;
	}
	if (fSPEX.DelaySigma <= 0) {
		PosDelaySigma	= false;
		SPEDelaySample	= fSPEX.Delay;
	}
	
	//Go along all arrival times
	for (int i = 0; i < int(ftimeseq->size()); i++) {		
		Interaction = true;
		
		if (fRND.Rndm() < fSPEX.Pphe_PC) {
		//If photon created any phe in PC
			AmplSigma = fSPEX.AmplSigma;
			
			if (fRND.Rndm() < fSPEX.DPE_PC) {
			//If photon created 2 phe in PC				
				Ampl = fSPEX.Ampl_DPE_PC;	//Amplitude from double photoelectron
				NPhe->Fill(2);	//Fill histogram with 2phe
			//	cout << i << " photon created 2 phe in PC" << endl;
			}
			else {
			//else photon created 1 phe in PC
				Ampl = fSPEX.Ampl;	//Single photoelectron amplitude
				NPhe->Fill(1);	//Fill histogram with 1phe
			//	cout << i << " photon created 1 phe in PC" << endl;
			}
		}
		else {
		//else check for interact with 1 dynode
			if (fRND.Rndm() < fSPEX.Pphe_1d) {
			//If photon created phe in 1dyn
				AmplSigma = fSPEX.AmplSigma_1d;
				Ampl = fSPEX.Ampl * fSPEX.ArbAmpl_1d;
				NPhe->Fill(3);	//"3" respond to 1phe from 1 dynode
			//	cout << i << " photon created 1 phe in 1 dynode " << endl;
			}
			else {
			//else photon didn't interact with PMT
				Interaction = false;
				NPhe->Fill(0);	//Fill histogram with no phe
			//	cout << i << " photon didn't create any phe" << endl;
			}
		}
		
		if (Interaction) {
			//If sigma for ampl or delay >0 , simulate
			if (PosAmplSigma)
				SPEAmplSample	= fRND.Gaus(Ampl, AmplSigma);
			if (PosDelaySigma)
				SPEDelaySample	= fRND.Gaus(fSPEX.Delay, fSPEX.DelaySigma);
			
			//For simulated amplitude fill PulseArea histogram
			PA->Fill(fSPE->GetParameter(2) * sqrt(2*TMath::Pi()) * abs(SPEAmplSample)/(ns*mV));
			
			//Calculate left and right samples including SPE
			arrtime = (*ftimeseq)[i] + SPEDelaySample - fOWX.Delay;
			startsample 	=  ceil((arrtime + fSPE->GetXmin())/fOWX.Period);
			finishsample	= floor((arrtime + fSPE->GetXmax())/fOWX.Period);
			if (startsample < 0)
				startsample = 0;
			if (finishsample > fOWX.Num - 1)
				finishsample = fOWX.Num - 1;			
		//	cout << "for " << i << " SPE AMPL=" << SPEAmplSample/mV << " mV";
		//	cout << " ARRTIME=" << (arrtime + fOWX.Delay)/ns << " ns" << endl;		
		//	cout << "Write from " << startsample << " to " << finishsample << " samples";
		//	cout << " (Domain = " << fSPEX.Domain << ")" << endl;
		
			//Add SPE to OutWave
			for (int sample = startsample; sample <= finishsample; sample++) {
				sampletime = sample*fOWX.Period;
				OutWave[sample] += SPEAmplSample*(fSPE->Eval(sampletime - arrtime))/fOWX.Gain;
			//	cout << "\tfor " << sample << " sample (" << sampletime/ns << " ns) added ";
			//	cout << SPEAmplSample*(fSPE->Eval(sampletime - arrtime))/fOWX.Gain << endl;
			}
		}
	}
//	cout << "OutWave was created" << endl;
}

//Old version for creating OutWave
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

//Print OutWave
void MakeWave::PrintOutWave() {
	for (int i = 0; i < fOWX.Num; i++) {
		cout << OutWave[i] << "\t(t = " << (fOWX.Delay + i*fOWX.Period)/ns << " ns)" << endl;
	}
	cout <<"OutWave was printed"<< endl;
}

//Draw OutWave
void MakeWave::Draw (TString name) {
	TCanvas *c1 = new TCanvas();
	TGraph *g1 = new TGraph(fOWX.Num);
	for (int i = 0; i < fOWX.Num; i++) {
		g1->SetPoint(i, (fOWX.Delay + i*fOWX.Period)/ns, OutWave[i]);
	}
	g1->SetTitle("Waveform");
	g1->Draw();
	if (name != "")
		c1->SaveAs(name);	//graph is closed when enabled
	c1->WaitPrimitive();
}