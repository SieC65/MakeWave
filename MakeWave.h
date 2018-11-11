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

class MakePhotons
{
	public:
	
		MakePhotons ();
		
		// SETTERS
		void SetTau (Double_t TauFast, Double_t TauSlow); // Set tau for fast & slow Sc exp decay functions
		void SetFastFrac (TF1* FastER_func, TF1* FastNR_func);  // Set functions of fast Sc fractions depend on photons number for ER and NR
		void SetFastFrac (Double_t FastER, Double_t FastNR); // Set constant fast Sc fractions
		void SetDefFastFract (); // Set default values for fast fraction functions
		
		// GETTERS
		vector <double> GetSimPhotonTimes() {return fSimPhotonTimes;}
		
		// ACTIONS
		void SimulatePhotons(Int_t NumPhotons, Double_t FracNR);

	private:
	
		// Auxiliary
		
		enum func_type {constant, function} fFast_type; // Show if fast Sc fraction depends on photons number
		enum InterType {ER, NR} fInterType; // Interaction type
		Int_t fMinPhotons;  // Minimum photons number for default function of fast Sc fraction
		Int_t fMaxPhotons;  // Maximum photons number (right edge of function)
		void SimPhTimes(Int_t NumPhotons, InterType recoil, TF1* FastProbFunc, TF1* ExpDecay);
		
		// Physics
		Double_t fTauFast;  // Tau for fast scintillation pdf
		Double_t fTauSlow;  // Tau for slow scintillation pdf
		Double_t fFastER;   // Fraction of fast component for scintillation from ER
		Double_t fFastNR;   // The same for NR
		TF1* fFastER_func;  // Fraction of fast component for Sc from ER depends on photons number
		TF1* fFastNR_func;  // The same for NR
		
		// Output
		vector <double> fSimPhotonTimes; // Output vector of photons arrival times
};

#endif // MakeWave_H
