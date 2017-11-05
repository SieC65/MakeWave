#ifndef MakeWave_H
#define MakeWave_H
#include <vector>
#include <Rtypes.h>
#include <TF1.h>
#include "SystemOfUnits.h"
#include <TRandom.h>
#include <TString.h>

using namespace std;
using namespace CLHEP;

class MakeWave
{
	public:
		MakeWave ();
		//SPE Parameters
		struct SPEPar {
			Bool_t Type;			//0 - gaussian , 1 - square pulse
			Double_t Trig;			//Coeff: Trigger for SPE = Trig*SPEAmpl
			Double_t Width;			//FWHM of SPE gaussian (or width of square pulse)
			Double_t Ampl;			//Amplitude of SPE gaussian
			Double_t AmplSigma;		//Sigma for SPE amplitude
			Double_t Delay;			//Delay of SPE signal
			Double_t DelaySigma;	//Sigma for SPE delay
			Double_t Domain;		//Domain width of SPE
		};
		//OutWave Parameters
		struct OutWavePar {
			Int_t Num ;			//Number of samples in OutWave
			Double_t Delay;		//Delay from "0" of abs.time to "0" sample of OutWave
			Double_t Period;	//Time between samples of OutWave
			Double_t Gain;		//Units of ADC
		};		
		void SetSPE (TF1 *SPE, SPEPar SPEX);	//Set SPE parameters
		void SetParams (OutWavePar OWX);		//Set OutWave parameters
		void SetTimeSeq (vector <double> *Tseq);	//Set vector of SPE's arrival times
		void CreateOutWave ();		//Create OutWave
		void CreateOutWaveOld();	//Create OutWave (old algorithm)
		void PrintOutWave ();		//Print all values of output signal
		void Draw (TString filename = "");		//Draw OutWave (on screen or to file)
		vector <double> OutWave;	//Output waveform
		void SetRand(Bool_t isrnd = true) {fisrnd = isrnd;}	//Randomize rnd. gen. for SPE params
		Bool_t GetIsrnd () {return fisrnd;}		//Return, was randomize set or wasn't
	private:
		SPEPar fSPEX;
		OutWavePar fOWX;
		vector <double> *ftimeseq;
		TF1 *fSPE;
		TRandom fRND;	//Random generator
		Bool_t fisrnd;
};

#endif // MakeWave_H
