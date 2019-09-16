#ifndef MakeWave_H
#define MakeWave_H

#include <vector>

#include <Rtypes.h>
#include "SystemOfUnits.h"

#include <REDFile/File.hh>
#include <REDEvent/Event.hh>

#include "PMT_R11410.hh"

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Class for generate output waveform from vector of photon times          //
// It allows to create REDFile compatible to REDOffline                    //
// Also it has tools for calculating F90 parameter (without dark current)  //
// Class that corresponds to electronic read-out from PMT                  //
// It allows to create output waveform from vector of photon times         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

using std::vector;
using CLHEP::ns;

class MakeWave
{
	public:

		MakeWave ();

	// SETTERS
		void SetPMT (RED::PMT* pmt); // Set PMT object
		void SetOutWave (Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay); // Set OutWave parameters
		void SetDefaults (); // Set default OutWave parameters
		void SetPhotonTimes (vector <double> *PhotonTimes); // Set vector of photon arrival times
		
	// GETTERS
		vector <double> GetOutWave () {return fOutWave;} // Output waveform vector
		// Get outWave parameters
		Double_t GetPeriod ()         {return fPeriod;}      // Time between samples of OutWave
		Double_t GetGain ()           {return fGain;}        // ADC resolution
		Double_t GetNumSamples ()     {return fNumSamples;}  // Number of samples in OutWave
		Double_t GetDelay ()          {return fDelay;}       // Delay from "0" of abs.time (related to photons times) to the left edge of OutWave

		// Tools for calculate F90
		Double_t GetFrac (Double_t FracWindow = 90*ns, Double_t TotalWindow = 0); // Fraction of light in the first FracWindow of pulse (default 90*ns)
		                                                                          // If needed, total pulse width can be limited by TotalWindow (leave 0 otherwise)
		Int_t GetNumPE () {return fPhotoElectrons->size();} // Get number of photoelectrons emitted last run

		// REDFile activities
		RED::OutputFile* GetCurrentFile () {return fOutFile;} // Return pointer to existing file
		
	// ACTIONS
		void CreateOutWave ();   // Create OutWave
		RED::OutputFile* GetNewFile (const char *filename); // Create new REDFile
		void AddToFile(); // Add current fOutWave to new event in REDFile
		void CloseFile(); // Close REDFile

	// OUTPUT
		void PrintOutWave ();    // Print OutWave (all times & amplitudes)
		void DrawHists ();       // Draw some histograms
		void DrawOutWave ();     // Draw OutWave
		void SaveOutWave (const char *filename = "OW.root", const char *title = "OutWave"); // Save OutWave to file

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
		TH1F* fPulseAreaHist; // Area under pulse SPE (DPE)

		//File
		RED::OutputFile *fOutFile;
		RED::Event *fEvent;
		RED::Waveform *fWaveform;
		RED::RunInfo *fRunInfo;
		Int_t fNumEv;
};

#endif // MakeWave_H
