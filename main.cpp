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
	Double_t Trig	 		= 0.01;	//Coeff: Trigger for SPE = Trig*SPEAmpl
	Double_t SPEWidth 		= 10*ns;	//FWHM of SPE gaussian
	Double_t SPEAmpl 		= (-0.001)*volt;	//Amplitude of SPE gaussian
	Double_t SPEAmplSigma	= 0.2*0.001*volt;	//Sigma for SPE amplitude
	Double_t SPEDelay		= 20*ns;	//Delay of SPE signal
	Double_t SPEDelaySigma	= 5*ns;		//Sigma for SPE delay
	Int_t Num 				= 20;		//Number of samples in OutWave
	Double_t Delay 			= 20*ns;	//Delay from "0" of abs.time to "0" sample of OutWave
	Double_t Period 		= 2*ns;		//Time between samples of OutWave
	Double_t Gain 			= 0.125*(volt/1000);	//Units of ADC
	
	//Let's create example of sequence of SPE's arrival times
	Double_t TseqArr[] = {/* -5.1 , -3.0 , 1 , 2 , */22 /*, 24 , 44 , 48 , 68 , 74 , 94 , 102 , 122 , 132 , 152 , 164 */};
	Int_t TseqArrSize = sizeof(TseqArr) / sizeof(TseqArr[0]);
	vector <double> Tseq (TseqArrSize);
	for (int i = 0; i < TseqArrSize; i++) {
		TseqArr[i] *=ns;	//set times in ns
		Tseq[i] = TseqArr[i];
	}
	
	TF1 *SPE = new TF1("SPE","((x > [0]) && (x < [1]))*[2]");	//Square pulse
	SPE->SetParameter(0, 0);			//Left end of the range
	SPE->SetParameter(1, SPEWidth);		//Right end of the range
	SPE->SetParameter(2, SPEAmpl);		//Amplitude of square pulse
/*	TF1 *SPE = new TF1("SPE","gaus(0)",0,(2*HalfDom(Trig, SPEWidth)));	//Gaussian
	SPE->SetParameter(0, SPEAmpl);		//Amplitude of gauss
	SPE->SetParameter(1, (HalfDom(Trig, SPEWidth)));	//Center
	SPE->SetParameter(2, (SPEWidth)/(2*sqrt(2*log(2))));	//Sigma
*/	
	MakeWave *a = new MakeWave;		//Create object of MakeWave class
	a->SetSPE (SPE, SPEAmpl, SPEAmplSigma, SPEDelay, SPEDelaySigma);	//Set SPE
	a->SetParams (Num, Delay, Period, Gain);	//Set parameters for OutWave
	a->SetTimeSeq (&Tseq);			//Set sequence of SPE arrival times
	a->SetIsrnd(true);
	cout << a->GetIsrnd() << endl;
	// Create average SPE
	vector <double> MeanOutWave;
	MeanOutWave.resize(Num);
	for (int i = 1; i < 500; i++) {
		a->CreateOutWave();				//Create OutWave vector
	//	a->Draw();
		for (int k = 0; k < Num; k++) {
			MeanOutWave[k] += a->OutWave[k];
		}
	}
	TCanvas *c2 = new TCanvas();
	TGraph *g2 = new TGraph(Num);
	for (int i = 0; i < Num; i++) {
		g2->SetPoint(i, (Delay + i*Period)/ns, MeanOutWave[i]);
	}
	g2->Draw();
	c2->WaitPrimitive();
		
	a->CreateOutWave();	//Create OutWave vector
	a->Draw();
	return 0;
}
