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

		// SETTERS
		void SetPMT (RED::PMT* PMT); // Set PMT object
		void SetOutWave (Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay); // Set OutWave parameters
		void SetPhotonTimes (vector <double> *PhotonTimes); // Set vector of photon arrival times
		

		// ACTIONS
		void CreateOutWave ();   // Create OutWave
		void test(RED::PMT_R11410 &); // Friend of PMT_R11410, draws hists

		// OUTPUT
		void PrintOutWave ();    // Print all values of output signal
		void DrawOutWave ();     // Draw OutWave

		// GETTERS
		vector <double> GetOutWave() {return fOutWave;} // Output waveform vector

	private:
	
		// OutWave Parameters
		vector <double> fOutWave;
		Double_t fPeriod;  // Time between samples of OutWave
		Double_t fGain;    // ADC resolution
		Int_t fNumSamples; // Number of samples in OutWave
		Double_t fDelay;   // Delay from "0" of abs.time (related to photons times) to "0" sample of OutWave
		
		// PMT
		RED::PMT *fPMT; // PMT object
		RED::PMT::PulseArray *PhotonPulse; // array of SPE caused by photons
		RED::PMT::PulseArray *DarkPulse;   // array of SPE caused by dark counts
		
		vector <double> *fPhotonTimes; // Photons arrival times
		void AddPulseArray (RED::PMT::PulseArray *Pulses); // Adding pulses to output waveform
};

#endif // MakeWave_H
