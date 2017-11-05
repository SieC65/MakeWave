#include <iostream>
#include "MakeWave.h"
#include <TF1.h>
#include <TApplication.h>

int main() {
	//Set parameters values
	TApplication *app = new TApplication("Canvas",0,0);
	TF1 *SPE = new TF1("SPE","((x > 0) && (x < 10))*(-0.001)");	//Function SPE, set in ns and mV
	Int_t Num 			= 10;			//Number of samples in OutWave
	Double_t Delay 		= 0 *ns;		//Delay from "0" of abs.time to "0" sample of OutWave
	Double_t Period 	= 2 *ns;		//Time between samples of OutWave
	Double_t Gain 		= 0.125*(volt/1000);	//Units of ADC
	Double_t SPEUnit 	= volt;			//Units of SPE
	//Let's create example of sequence of SPE's arrival times
		Double_t TseqArr[] = { -5.1 , -3.0 , 1.6 , 1.8 , 1.9 , 4.3 , 7.1 , 7.5 };
		Int_t TseqArrSize = sizeof(TseqArr) / sizeof(TseqArr[0]);
		vector <double> Tseq (TseqArrSize);
		for (int i = 0; i < TseqArrSize; i++) {
			TseqArr[i] *=ns;	//set times in ns
			Tseq[i] = TseqArr[i];
		}
	MakeWave *a = new MakeWave;			//Create object of MakeWave class
	a->SetSPE (SPE);					//Set SPE function
	a->SetParams (Num, Delay, Period, Gain, SPEUnit);	//Set parameters
	a->SetTimeSeq (&Tseq);				//Set sequence of SPE's arrival times
	a->CreateOutWave();					//Create OutWave vector
	a->PrintOutWave();					//Print OutWave
	a->Draw();
	return 0;
}
