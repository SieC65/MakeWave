#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TApplication.h>
#include <TMath.h>

#include <TCanvas.h>
#include <TGraph.h>

//Time between center of gaus and trigger moment. Domain = 2*HalfDom.
Double_t HalfDom (Double_t k, Double_t Width) {
	Double_t HaD = Width*sqrt(log2(1/k))/2;
	return HaD;
}

int main() {
	TApplication *app = new TApplication("Canvas",0,0);
	
	//Set parameters values
	Bool_t SPEType 			= 1;		//0 - gaussian , 1 - square pulse
	Double_t Trig	 		= 0.01;		//Coeff: Trigger for SPE = Trig*SPEAmpl
	Double_t SPEWidth 		= 10*ns;	//FWHM of SPE gaussian (or width of square pulse)
	Double_t SPEAmpl 		= (-0.001)*volt;	//Amplitude of SPE gaussian
	Double_t SPEAmplSigma	= 0.2*0.001*volt;	//Sigma for SPE amplitude
	Double_t SPEDelay		= 20*ns;	//Delay of SPE signal
	Double_t SPEDelaySigma	= 0.05*ns;	//Sigma for SPE delay
	Int_t Num 				= 50;		//Number of samples in OutWave
	Double_t Delay 			= 0*ns;		//Delay from "0" of abs.time to "0" sample of OutWave
	Double_t Period 		= 2*ns;		//Time between samples of OutWave
	Double_t Gain 			= 0.125*(volt/1000);	//Units of ADC
	
	//Set sequence of SPE arrival times
	Double_t TseqArr[] = { -5.1 , -3.0 , 1 , 2 , 24 , 44 , 48 , 68 , 74 , 94 , 102 , 122 , 132 , 152 , 164 };
	Int_t TseqArrSize = sizeof(TseqArr) / sizeof(TseqArr[0]);
	vector <double> Tseq (TseqArrSize);
	for (int i = 0; i < TseqArrSize; i++) {
		TseqArr[i] *=ns;	//Set times in ns
		Tseq[i] = TseqArr[i];
	}
	
	//SPE form is set below
	Double_t SPEDomain;	//Domain width of SPE
	TF1 *SPE;
	if (SPEType == 1) {
		//Square pulse
		SPEDomain = SPEWidth;
		SPE = new TF1("SPE","((x > [0]) && (x < [1]))*[2]");
		SPE->SetParameter(0, 0);			//Left end of the range
		SPE->SetParameter(1, SPEDomain);		//Right end of the range
		SPE->SetParameter(2, SPEAmpl);		//Amplitude of square pulse
	}
	else {
		//Gaussian
		SPEDomain = 2*HalfDom(Trig, SPEWidth);
		SPE = new TF1("SPE","gaus(0)",0,(2*HalfDom(Trig, SPEWidth)));
		SPE->SetParameter(0, SPEAmpl);		//Amplitude of gauss
		SPE->SetParameter(1, SPEDomain/2);	//Center
		SPE->SetParameter(2, (SPEWidth)/(2*sqrt(2*log(2))));	//Sigma
	}

	//Set parameters given above and create OutWave
	MakeWave *a = new MakeWave;		//Create object of MakeWave class
	a->SetSPE (SPE, SPEAmpl, SPEAmplSigma, SPEDelay, SPEDelaySigma, SPEDomain);	//Set SPE parameters
	a->SetParams (Num, Delay, Period, Gain);	//Set OutWave parameters
	a->SetTimeSeq (&Tseq);			//Set sequence of SPE arrival times
	a->SetIsrnd(true);				//Randomize rnd. gen. for SPE params
	a->CreateOutWave();				//Create OutWave vector
//	a->Draw("one.jpg");				//Draw OutWave on screen and to file

	
	//DEBUG BEGIN
	//Create average OutWave for NMAX OutWaves with the same SPE arrival times
	Int_t NMAX = 200;
	vector <double> MeanOutWave;
	vector <double> MeanOutWaveOld;
	MeanOutWave.resize(Num);
	MeanOutWaveOld.resize(Num);
	for (int i = 1; i < NMAX; i++) {
		a->CreateOutWave();
		for (int k = 0; k < Num; k++) {
			MeanOutWave[k] += a->OutWave[k];
		}
		a->CreateOutWaveOld();
		for (int k = 0; k < Num; k++) {
			MeanOutWaveOld[k] += a->OutWave[k];
		}
	}
	TCanvas *c2 = new TCanvas();
	c2->Divide(1,2);
	
	c2->cd(1);
	TGraph *g2 = new TGraph(Num);
	for (int i = 0; i < Num; i++) {
		g2->SetPoint(i, (Delay + i*Period)/ns, MeanOutWave[i]);
	}
	g2->SetTitle("New");
	g2->Draw();
	
	c2->cd(2);
	TGraph *g3 = new TGraph(Num);
	for (int i = 0; i < Num; i++) {
		g3->SetPoint(i, (Delay + i*Period)/ns, MeanOutWaveOld[i]);
	}
	g3->SetTitle("Old");
	g3->Draw();
	
//	c2->SaveAs("ave.jpg");	//graph is closed if uncommented
	c2->WaitPrimitive();
	//DEBUG END
	
	return 0;
}
