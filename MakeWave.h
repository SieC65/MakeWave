#ifndef MakeWave_H
#define MakeWave_H
#include <vector>
#include <Rtypes.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TString.h>
#include "SystemOfUnits.h"
#include <TH1.h>

namespace CLHEP {
	static const double mV = 1.e-3*volt;
}

using namespace std;

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
			Double_t QE;			//Quantum Efficiency (full)
			Double_t DPE_PC;		//Double Photoelectron Emission = P(2phe)/(P(2phe)+P(1phe)) for PC
			Double_t QE_1d;			//Quantum Efficiency for 1dyn (only 1phe)
			Double_t ArbCurr_1d;	//Current(1d)/Current(PC) for SPE
			Double_t ArbAmpl_1d;	//Ampl(1d)/Ampl(PC) for SPE
			Double_t AmplSigma_1d;	//Sigma for amplitude from 1dyn
			Double_t GF_1d;			//Arb. geometrical factor of 1dyn
			Double_t QE_1d_ratio;	//auxiliary quantity
			Double_t Pphe_PC;		//Probability of photoeffect on PC
			Double_t Pphe_1d;		//Probability of photoeffect on 1dyn for photons passed PC
			Double_t Ampl_DPE_PC;	//Amplitude of signal from 2phe from PC
		};
		//OutWave Parameters
		struct OutWavePar {
			Double_t Period;	//Time between samples of OutWave
			Double_t Gain;		//Units of ADC
			Int_t Num ;			//Number of samples in OutWave
			Double_t Delay;		//Delay from "0" of abs.time to "0" sample of OutWave
		};		
		void SetSPE (const SPEPar &SPEX);	//Set SPE parameters
		void SetParams (const OutWavePar &OWX);	//Set OutWave parameters
		void SetTimeSeq (vector <double> *Tseq);	//Set vector of photon arrival times
		void CreateOutWave ();		//Create OutWave
		void CreateOutWaveOld();	//Create OutWave (old algorithm)
		void PrintOutWave ();		//Print all values of output signal
		void Draw (TString filename = "");	//Draw OutWave (on screen or to file <filename>)
		vector <double> OutWave;	//Output waveform
		void SetRand(Bool_t isrnd = true);	//Randomize rnd. gen. for SPE params
		Bool_t GetIsrnd ();					//Return, was randomize set or wasn't
		TH1F *NPhe;					//Number of photons created 0, 1 and 2 phe
		TH1F *PA;					//Area under pulse SPE (DPE)
	private:
		SPEPar fSPEX;
		OutWavePar fOWX;
		vector <double> *ftimeseq;
		TF1 *fSPE;
		TRandom3 fRND;	//Random generator
		Bool_t fisrnd;
};

#endif // MakeWave_H
