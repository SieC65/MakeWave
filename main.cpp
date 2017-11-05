#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TApplication.h>
#include <TMath.h>

//Time between center of gaus and trigger moment. Domain = 2*HalfDom.
Double_t HalfDom (Double_t k, Double_t Width) {
	Double_t HaD = Width*sqrt(log2(1/k))/2;
	return HaD;
}

int main() {
	TApplication *app = new TApplication("Canvas",0,0);
	//Set parameters values
	Double_t trig	 		= 0.2;		//trigger for SPE = trig*SPEAmplitude
	Double_t SPEamplitude 	= (-0.001)*volt;	//Amplitude of SPE gaussian
	Double_t SPEwidth 		= 10*ns;	//FWHM of SPE gaussian
	Int_t Num 				= 100;		//Number of samples in OutWave
	Double_t Delay 			= 0*ns;		//Delay from "0" of abs.time to "0" sample of OutWave
	Double_t Period 		= 2*ns;		//Time between samples of OutWave
	Double_t Gain 			= 0.125*(volt/1000);	//Units of ADC
	//Let's create example of sequence of SPE's arrival times
	Double_t TseqArr[] = { -5.1 , -3.0 , 1 , 2 , 22 , 24 , 44 , 48 , 68 , 74 , 94 , 102 , 122 , 132 , 152 , 164};
		Int_t TseqArrSize = sizeof(TseqArr) / sizeof(TseqArr[0]);
		vector <double> Tseq (TseqArrSize);
		for (int i = 0; i < TseqArrSize; i++) {
			TseqArr[i] *=ns;	//set times in ns
			Tseq[i] = TseqArr[i];
		}
	TF1 *SPE = new TF1("SPE","((x > [0]) && (x < [1]))*[2]");	//Square pulse
	SPE->SetParameter(0, 0);				//Left end of the range
	SPE->SetParameter(1, SPEwidth);			//Right end of the range
	SPE->SetParameter(2, SPEamplitude);	//Amplitude of square pulse
/*	TF1 *SPE = new TF1("SPE","gaus(0)",0,(2*HalfDom(trig, SPEwidth)));	//Gaussian
	SPE->SetParameter(0, SPEamplitude);		//Amplitude of gauss
	SPE->SetParameter(1, (HalfDom(trig, SPEwidth)));	//Center
	SPE->SetParameter(2, (SPEwidth)/(2*sqrt(2*log(2))));	//Sigma
*/	SPE->Draw();
	MakeWave *a = new MakeWave;			//Create object of MakeWave class
	a->SetSPE (SPE);					//Set SPE function
	a->SetParams (Num, Delay, Period, Gain);	//Set parameters
	a->SetTimeSeq (&Tseq);				//Set sequence of SPE's arrival times
	a->CreateOutWave();					//Create OutWave vector
//	a->PrintOutWave();					//Print OutWave
	a->Draw();
	return 0;
}
