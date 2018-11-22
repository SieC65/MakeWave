#ifndef MakeWave_H
#define MakeWave_H

#include <vector>

#include <Rtypes.h>
#include <TString.h>

#include "PMT_R11410.hh"
#include <REDFile/File.hh>
#include <REDEvent/Event.hh>

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
		void SetDefaults ();
		

		// ACTIONS
		void CreateOutWave ();   // Create OutWave
		void DrawHists();        // Draw histograms

		// OUTPUT
		void PrintOutWave ();    // Print all values of output signal
		void DrawOutWave ();     // Draw OutWave
		void SaveOutWave (const char *filename); // Save OutWave to REDOffline-compatible *.root file

		// GETTERS
		vector <double> GetOutWave() {return fOutWave;} // Output waveform vector
		Double_t GetPeriod()         {return fPeriod;}
		Double_t GetGain()           {return fGain;}
		Double_t GetNumSamples()     {return fNumSamples;}
		Double_t GetDelay()          {return fDelay;}

	private:
	
		// OutWave Parameters
		vector <double> fOutWave;
		Double_t fPeriod;     // Time between samples of OutWave
		Double_t fGain;       // ADC resolution
		Int_t    fNumSamples; // Number of samples in OutWave
		Double_t fDelay;      // Delay from "0" of abs.time (related to photons times) to "0" sample of OutWave
		
		// PMT
		RED::PMT *fPMT; // PMT object
		RED::PMT::PulseArray *fPhotoElectrons;  // array of SPE caused by photons
		RED::PMT::PulseArray *fDarkElectrons;   // array of SPE caused by dark counts
		
		vector <double> *fPhotonTimes; // Photons arrival times
		void AddPulseArray (RED::PMT::PulseArray *Pulses); // Adding pulses to output waveform

		// Histograms
		TH1F* fDarkTimeHist;  // Time of dark pulses
		TH1F* fDarkAmplHist;  // Amplitude of dark pulses
		TH1F* fNumPheHist;    // Number of photoelectron created by 1 photon:
							 // -1 or -2: 1 or 2 phe in 1dyn
							 // +1 or +2: 1 or 2 phe in PC
		TH1F* fPulseAreaHist; // Area under pulse SPE (DPE)
};

#endif // MakeWave_H
