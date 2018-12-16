#ifndef SimPhotons_H
#define SimPhotons_H

#include <iostream>
#include <algorithm>
#include <vector>

#include <Rtypes.h>

#include <TH1.h>
#include <TF1.h>
#include <TSpline.h>
#include <TRandom3.h>
#include "SystemOfUnits.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "MakeWave.h"

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
// Class for simulate time moments of flashing photons due to 1 interaction event.   //
//                                                                                   //
// User should set parameters of scintillation caused by interaction                 //
// either as two numbers of photons corresponding to fast & slow scintillation       //
// or as a total number of emitted photons + type of interaction ("ER" | "NR").      //
//                                                                                   //
// After setting parameters of scintillation and interaction                         //
// user can get vector of photons times, where interaction happened at "0" point.    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

using std::vector;
using CLHEP::ns;

class SimPhotons
{
	public:
	
		SimPhotons ();
		
	// SETTERS
		void SetTau (Double_t TauFast, Double_t TauSlow); // Set tau for fast & slow scintillation exp(-t/tau) functions

		// Set fast fractions of scintillation for ER and NR processes
		void SetFastFrac (TF1* FastER_func, TF1* FastNR_func);  // as functions depend on photons number
		void SetFastFrac (Double_t FastER, Double_t FastNR);    // as constants
		void SetDefFastFract (); // Set fast fractions as (A + B / NumPhotons) with default A,B for ER & NR
		
	// GETTERS
		vector <double> GetSimPhotonTimes() {return fSimPhotonTimes;} // return vector of photons times
		
	// ACTIONS
		// Simulate flashing times
		vector <double> SimulatePhotons (Int_t NumPhotons, Double_t FastFrac);
		vector <double> SimulatePhotons (Int_t NumFast, Int_t NumSlow);
		vector <double> SimulatePhotons (Int_t NumPhotons, Option_t* type);

	private:
	
		// Auxiliary		
		enum func_type {constant, function} fFast_type; // Show if fast Sc fraction depends on photons number
		enum InterType {ER, NR} fInterType; // Interaction type
		Int_t fMinPhotons;  // Minimum photons number for default function of fast Sc fraction
		Int_t fMaxPhotons;  // Maximum photons number (right edge of function)
		
		// Physics
		Double_t fTauFast;  // Tau for fast scintillation pdf
		Double_t fTauSlow;  // Tau for slow scintillation pdf
		Double_t fFastER;   // Fraction of fast component for scintillation from ER
		Double_t fFastNR;   // The same for NR
		TF1* fFastER_func;  // Fraction of fast component for Sc from ER depends on photons number
		TF1* fFastNR_func;  // The same for NR
		TF1* fExpDecaySlow; // decay PDF of slow Sc component
		TF1* fExpDecayFast; // decay PDF of fast Sc component
		TF1* fFastProbFunc; // binomial PDF for fast fraction
		
		// Output
		vector <double> fSimPhotonTimes; // Output vector of photons arrival times
};

#endif // SimPhotons_H 
