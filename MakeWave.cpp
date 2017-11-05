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

//Set PMT parameters
void MakeWave::SetPMT (PMT *PMTX) {
	fPMT = PMTX;
}

//Set OutWave parameters
void MakeWave::SetOutWave (const OutWavePar &OWX) {
	fOWX = OWX;
	cout << "Parameters were set" << endl;
}

//Set sequence of photon times
void MakeWave::SetTimeSeq (vector <double> *Tseq) {
	ftimeseq = Tseq;	//Now addresses of ftimeseq and Tseq are equal (pointer Tseq created in main.cpp)
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
	Double_t sampletime = 0;	//Time of sample from "0" of OutWave
	Double_t arrtime	= 0;	//Arrival time from "0" of OutWave
	OutWave.clear();				//Clear vector OutWave
	OutWave.resize (fOWX.Num, 0);	//Clear vector OutWave
	Double_t SPEAmplSample	= 0;	//Random SPE amplitude
	Double_t SPEDelaySample	= 0;	//Random SPE delay
	Int_t startsample		= 0;	//First sample in SPE domain
	Int_t finishsample		= 0;	// Last sample in SPE domain
	Double_t Ampl			= 0;	//Amplitude of signal from photon
	Double_t AmplSigma		= 0;	//Sigma of amplitude
	Double_t Delay			= 0;	//Delay of SPE
	Int_t NumPhe			= 0;	//Number of emitted photoelectrons
	Double_t PA_Sum 		= 0;	//PulseArea for photon interaction
	Int_t InteractType		= 0;	//Interaction Type: 0 - no interaction;
									//1 or 2 - for Nphe in PC, 3 - for 1phe in 1dyn
		
	//Go along all arrival times
	for (int i = 0; i < int(ftimeseq->size()); i++) {		
		ShootPhoton (&NumPhe, &Ampl, &AmplSigma, &InteractType, &Delay);	//Define interaction type
		PA_Sum = 0;	//Initialize PulseArea value
		
		for (int k=0; k<NumPhe; k++) {
			
			//If sigma for ampl or delay >0 , spread it
			if (fSPEX.AmplSigma > 0)
				SPEAmplSample	= fRND.Gaus(Ampl, AmplSigma);
			else
				SPEAmplSample	= Ampl;
			if (fSPEX.DelaySigma > 0)
				SPEDelaySample	= fRND.Gaus(Delay, fSPEX.DelaySigma);
			else
				SPEDelaySample	= Delay;
			
			//Set value for add it to PulseArea histogram
			//In case of 2phe this value will first increase before addition
			PA_Sum += fSPE->GetParameter(2) * sqrt(2*TMath::Pi()) * abs(SPEAmplSample)/(ns*mV);

			//Calculate left and right samples including SPE
			arrtime = (*ftimeseq)[i] + SPEDelaySample - fOWX.Delay;
			startsample 	=  ceil((arrtime + fSPE->GetXmin())/fOWX.Period);
			finishsample	= floor((arrtime + fSPE->GetXmax())/fOWX.Period);
			if (startsample < 0)
				startsample = 0;
			if (finishsample > fOWX.Num - 1)
				finishsample = fOWX.Num - 1;
			
		/*	if (InteractType==3) {
				cout << "It was interaction in 1st dynode at " << (arrtime + fOWX.Delay)/ns << " ns" << endl;
				cout << "Ampl= ";
			}
		*/
			//Add SPE to OutWave
			for (int sample = startsample; sample <= finishsample; sample++) {
				sampletime = sample*fOWX.Period;
				OutWave[sample] += SPEAmplSample*(fSPE->Eval(sampletime - arrtime))/fOWX.Gain;
			}
		}
		
		if 	(!!NumPhe) {
			//For simulated amplitude fill PulseArea histogram
			PA->Fill(PA_Sum);
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
		if (fisrnd == true) {
			fRND.SetSeed(0);
		}
		SPEAmplVec[i] 	= fRND.Gaus(1, fSPEX.AmplSigma/fSPEX.Ampl);
		SPEDelayVec[i] 	= fRND.Gaus(fSPEX.Delay, fSPEX.DelaySigma);
	}
	for (int sample = 0; sample < fOWX.Num; sample++) {
		t = fOWX.Delay + sample*fOWX.Period;
		for (int i = 0; i < int(ftimeseq->size()); i++) {
			OutWave[sample] += (SPEAmplVec[i])*(fSPE->Eval(t - SPEDelayVec[i] - (*ftimeseq)[i]))/fOWX.Gain;
		}
	}
//	cout << "OutWave was created"<< endl;
}

//Define interaction type and, respectively, set Ampl, Ampl. sigma and Number of phe
void MakeWave::ShootPhoton (Int_t *NumPhe, Double_t *Ampl, Double_t *AmplSigma, Int_t *InteractType, Double_t *Delay) {
	Double_t RND	= 0;	//Random in range 0..1 defining interact. type
	RND = fRND.Rndm();
	
	if (RND < fPMT->GetPphe_PC()) {
	//If photon created any phe in PC
		*AmplSigma	= fSPEX.AmplSigma;
		*Ampl		= fSPEX.Ampl;
		*Delay		= fSPEX.Delay;
		
		if (RND < fPMT->GetDPE() * fPMT->GetPphe_PC()) {
		//If photon created 2 phe in PC				
			*NumPhe			= 2;
			*InteractType	= 2;
		}
		else {
		//else photon created 1 phe in PC
			*NumPhe			= 1;
			*InteractType	= 1;
		}
	}
	else {
	//else check for interact with 1 dynode
		if (RND > 1 - (1 - fPMT->GetPphe_PC()) * fPMT->GetPphe_1d()) {
		//If photon created phe in 1dyn
			*AmplSigma = fSPEX.AmplSigma * fPMT->GetArbAmpl_1d();
			*Ampl	= fSPEX.Ampl * fPMT->GetArbAmpl_1d();
			*Delay	= fSPEX.Delay - fPMT->GetTOF_1d();
			*NumPhe			= 1;
			*InteractType	= 3;
		}
		else {
		//else photon didn't interact with PMT
			*NumPhe			= 0;
			*InteractType	= 0;
		}
	}
	NPhe->Fill(*InteractType);	//Fill histogram with Nphe
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
		c1->SaveAs(name);	//Graph is closed when enabled
//	c1->WaitPrimitive();	
}