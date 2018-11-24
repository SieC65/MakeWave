#ifndef MakeWave_H
#define MakeWave_H

#include <vector>

#include <Rtypes.h>

#include "PMT_R11410.hh"

// Class that corresponds to electronic read-out from PMT
// It allows to create output waveform from vector of photon times

using std::vector;

class MakeWave
{
	public:

		MakeWave ();

		// SETTERS
		void SetPMT (RED::PMT* pmt); // Set PMT object
		void SetOutWave (Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay); // Set OutWave parameters
		void SetDefaults (); // Set default OutWave parameters
		void SetPhotonTimes (vector <double> *PhotonTimes); // Set vector of photon arrival times
		
		// ACTIONS
		void CreateOutWave ();   // Create OutWave

		// OUTPUT
		void PrintOutWave ();    // Print OutWave (all times & amplitudes)
		void DrawHists ();       // Draw some histograms
		void DrawOutWave ();     // Draw OutWave
		void SaveOutWave (const char *filename); // Save OutWave to REDOffline-compatible *.root file

		// GETTERS
		vector <double> GetOutWave () {return fOutWave;} // Output waveform vector
		// OutWave parameters
		Double_t GetPeriod ()         {return fPeriod;}      // Time between samples of OutWave
		Double_t GetGain ()           {return fGain;}        // ADC resolution
		Double_t GetNumSamples ()     {return fNumSamples;}  // Number of samples in OutWave
		Double_t GetDelay ()          {return fDelay;}       // Delay from "0" of abs.time (related to photons times) to the left edge of OutWave

	private:
		
		// FUNCTIONS
		void AddPulseArray (RED::PMT::PulseArray *Pulses); // Adding pulses to output waveform

		// VALUES

		// OutWave Parameters
		vector <double> fOutWave;
		Double_t fPeriod;
		Double_t fGain;
		Int_t    fNumSamples;
		Double_t fDelay;
		
		// PMT
		RED::PMT *fPMT; // PMT object
		RED::PMT::PulseArray *fPhotoElectrons;  // array of SPE caused by photons
		RED::PMT::PulseArray *fDarkElectrons;   // array of SPE caused by dark counts
		
		vector <double> *fPhotonTimes; // Photons arrival times

		// Histograms
		TH1F* fDarkTimeHist;  // Time of dark pulses
		TH1F* fDarkAmplHist;  // Amplitude of dark pulses
		TH1F* fNumPheHist;    // Number of photoelectron created by 1 photon:
							  // -1 or -2: 1 or 2 phe in 1dyn
							  // +1 or +2: 1 or 2 phe in PC
		TH1F* fPulseAreaHist; // Area under pulse SPE (DPE)
};

#endif // MakeWave_H
