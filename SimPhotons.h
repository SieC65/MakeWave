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

//using namespace RED;
using std::vector;
using CLHEP::ns;

class SimPhotons
{
	public:
	
		SimPhotons ();
		
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

#endif // SimPhotons_H 
