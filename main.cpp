#include <iostream>
#include "MakeWave.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "SystemOfUnits.h"

using namespace CLHEP;

//Print 2 graphs (new and old algorithms) for sum of DebugN OutWaves
void Compare (int DebugN, struct MakeWave::OutWavePar OW, MakeWave *aex);

int main() {
	TApplication *app = new TApplication("Canvas",0,0);
	Int_t DebugAlg	= 1;	//0 - for only new algorithm, any other natural - for Compare ()
	
	//Set SPE and OutWave parameters
	MakeWave::SPEPar SPEEx;		//SPE Example
		SPEEx.Type 			= 0;				//0 - gaussian , 1 - square pulse
		SPEEx.Trig	 		= 0.01;				//Coeff: Trigger for SPE = Trig*SPEAmpl
		SPEEx.Width 		= 10*ns;			//FWHM of SPE gaussian (or width of square pulse)
		SPEEx.Ampl 			= (-0.001)*volt;	//Amplitude of SPE gaussian
		SPEEx.AmplSigma		= 0.2*0.001*volt;	//Sigma for SPE amplitude
		SPEEx.Delay			= 20*ns;			//Delay of SPE signal
		SPEEx.DelaySigma	= 5*ns;				//Sigma for SPE delay
		SPEEx.Domain		= SPEEx.Width;		//Domain width of SPE
	MakeWave::OutWavePar OWEx;	//OutWave example
		OWEx.Num 		= 500;					//Number of samples in OutWave
		OWEx.Delay 		= 0*ns;					//Delay from "0" of abs.time to "0" sample of OutWave
		OWEx.Period 	= 0.2*ns;					//Time between samples of OutWave
		OWEx.Gain 		= 0.125*(volt/1000);	//Units of ADC
	
	//Set sequence of SPE arrival times
	Double_t TseqArr[] = { -20 , 0 , 5.5 , 15 , 50};
	Int_t TseqArrSize = sizeof(TseqArr) / sizeof(TseqArr[0]);
	vector <double> Tseq (TseqArrSize);
	for (int i = 0; i < TseqArrSize; i++) {
		TseqArr[i] *=ns;	//Set times in ns
		Tseq[i] = TseqArr[i];
	}

	//Set parameters given above and create OutWave
	MakeWave *a = new MakeWave;	//Create object of MakeWave class
	a->SetSPE (SPEEx);			//Set SPE
	a->SetParams (OWEx);		//Set OutWave parameters
	a->SetTimeSeq (&Tseq);		//Set sequence of SPE arrival times
	a->SetRand(true);			//Randomize random generator for SPE params
	if (DebugAlg == 0) {
		a->CreateOutWave();		//Just create OutWave vector
		a->Draw();				//And draw OutWave on screen (may to file)
	}
	else {
		Compare (DebugAlg, OWEx, a);	//Debug algorithm
	}
	return 0;
}

//DEBUG FUNC
//Create average OutWave for DebugN OutWaves with the same SPE arrival times
void Compare (int DebugN, struct MakeWave::OutWavePar OW, MakeWave *aex) {
	vector <double> MeanOutWaveNew;		//sum of DebugN OutWaves from new alg.
	vector <double> MeanOutWaveOld;		//sum of DebugN OutWaves from old alg.
	MeanOutWaveNew.resize(OW.Num);
	MeanOutWaveOld.resize(OW.Num);
	for (int i = 1; i <= DebugN; i++) {
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
//	c2->SaveAs("ave.jpg");	//graph is closed if uncommented
	c2->WaitPrimitive();
}