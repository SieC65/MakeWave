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
		//Set SPE parameters
		void SetSPE (TF1 *SPE, Double_t SPEAmpl, Double_t SPEAmplSigma, Double_t SPEDelay, Double_t SPEDelaySigma, Double_t SPEDomain);
		//Set OutWave parameters
		void SetParams (Int_t Num, Double_t Delay, Double_t Period, Double_t Gain);
		//Set vector of SPE's arrival times
		void SetTimeSeq (vector <double> *Tseq);
		void CreateOutWave ();		//Create OutWave
		void CreateOutWaveOld();	//Create OutWave (old)
		void PrintOutWave ();		//Print all values of output signal
		void Draw (TString filename = "");			//Draw OutWave values from abs time
		vector <double> OutWave;	//Output waveform
		void SetIsrnd(Bool_t isrnd = true) {fisrnd = isrnd;}	//Randomize rnd. gen. for SPE params
		Bool_t GetIsrnd () {return fisrnd;}		//Return, was randomize set or wasn't
	private:
		vector <double> *ftimeseq;
		Double_t fPeriod;
		Double_t fDelay;
		Double_t fGain;
		Int_t fNum;
		TF1 *fSPE;
		TRandom fRND;	//Random generator
		Bool_t fisrnd;
		Double_t fSPEAmpl;
		Double_t fSPEAmplSigma;
		Double_t fSPEDelay;
		Double_t fSPEDelaySigma;
		Double_t fSPEDomain;
};

#endif // MakeWave_H
