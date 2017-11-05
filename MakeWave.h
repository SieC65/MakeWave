#ifndef MakeWave_H
#define MakeWave_H

#include <vector>

#include <Rtypes.h>

#include "PMT_R11410.hh"

//using namespace RED;
using std::vector;

class MakeWave
{
	public:
		MakeWave ();
		// OutWave Parameters
		struct OutWavePar {
			Double_t Period; // Time between samples of OutWave
			Double_t Gain;   // Units of ADC
			Int_t Num ;      // Number of samples in OutWave
			Double_t Delay;  // Delay from "0" of abs.time to "0" sample of OutWave
		};      
		void SetPMT (RED::PMT* PMTX);
		void SetOutWave (const OutWavePar &OWX); // Set OutWave parameters
		void SetTimeSeq (vector <double> *Tseq); // Set vector of photon arrival times
		void CreateOutWave ();   // Create OutWave
		void PrintOutWave ();    // Print all values of output signal
		void Draw ();            // Draw OutWave
		vector <double> OutWave; // Output waveform
	private:
		RED::PMT *fPMT;
		OutWavePar fOWX;
		vector <double> *ftimeseq;
		RED::PMT::PulseArray *PhotonPulse;
		RED::PMT::PulseArray *DarkPulse;
		void AddPulseArray (RED::PMT::PulseArray *Pulses);
};

#endif // MakeWave_H
