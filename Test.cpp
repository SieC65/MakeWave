#include <iostream>

#include <TGraph.h>
#include "SystemOfUnits.h"

#include "Test.h"

using std::cout;
using std::endl;
using std::vector;
using CLHEP::ns;

MakePhotons::MakePhotons() {
	fMinPhotons = 140;
	fMaxPhotons = 2000;
	fTauFast    = 7*ns;
	fTauSlow    = 1600*ns;
	SetFastFrac(0.22, 0.75);
}

void MakePhotons::SetTau (Double_t TauFast, Double_t TauSlow) {
	fTauFast = TauFast;
	fTauFast = TauSlow;
}

void MakePhotons::SetFastFrac (TF1* FastER_func, TF1* FastNR_func) {
	fFast_type = function;
	fFastER_func = FastER_func;
	fFastNR_func = FastNR_func;
}

void MakePhotons::SetFastFrac (Double_t FastER, Double_t FastNR) {
	fFast_type = constant;
	fFastER = FastER;
	fFastNR = FastNR;	
}

void MakePhotons::SetDefFastFract () {
	fFast_type = function;
	fFastER_func = new TF1 ("Fast Sc Fraction for ER","[0]+[1]/x",fMinPhotons,fMaxPhotons);
	fFastER_func -> SetParameter (0, 0.178464);
    fFastER_func -> SetParameter (1, 46.705);
	fFastNR_func = new TF1 ("Fast Sc Fraction for NR","[0]+[1]/x",fMinPhotons,fMaxPhotons);
	fFastNR_func -> SetParameter (0, 0.723801);
    fFastNR_func -> SetParameter (1, -52.0528);
}

void MakePhotons::SimPhTimes(Int_t NumPhotons, InterType recoil, TF1* FastProbFunc, TF1* ExpDecay) {
	cout << "recoil type:" << recoil << endl;
	// Calculate numbers of fast & slow photons for each group
	Double_t FastProb = 0;
	if      ((fFast_type == function) && (recoil == ER))
	{
		FastProb = fFastER_func -> Eval(NumPhotons);
		cout << "ER case, function" << endl;
	}
	else if ((fFast_type == function) && (recoil == NR))
	{
		FastProb = fFastNR_func -> Eval(NumPhotons);
		cout << "NR case, function" << endl;
	}
	else if ((fFast_type == constant) && (recoil == ER))
	{
		FastProb = fFastER;
		cout << "ER case, const" << endl;
	}
	else if ((fFast_type == constant) && (recoil == NR))
	{
		FastProb = fFastNR;
		cout << "NR case, const" << endl;
	}
	FastProbFunc -> SetRange (0, NumPhotons+1);
	FastProbFunc -> SetParameter (0, FastProb);
	FastProbFunc -> SetParameter (1, NumPhotons);
	cout << recoil << ": fast fraction = " << FastProb << endl;
	Int_t NumFast = FastProbFunc -> GetRandom (0, NumPhotons+1);
	Int_t NumSlow = NumPhotons - NumFast;
	
	// Simulate their flashing times
	ExpDecay -> SetRange (0, 500*ns);
	ExpDecay -> SetParameter (0, fTauFast);
	for (Int_t phe = 0; phe < NumFast; phe++)
	{
		fSimPhotonTimes.push_back(ExpDecay->GetRandom(0,500*ns));
	}
	ExpDecay -> SetRange (0, 50000*ns);
	ExpDecay -> SetParameter (0, fTauSlow);
	for (Int_t phe = 0; phe < NumSlow; phe++) 
	{
		fSimPhotonTimes.push_back(ExpDecay->GetRandom(0,50000*ns));
	}
}

void MakePhotons::SimulatePhotons(Int_t NumPhotons, Double_t FracNR) {
	cout << "Start simulating photons" << endl;
	
	// First, divide photons to NR and ER groups
	Int_t NumPhotonsNR = round (NumPhotons*FracNR); // Number of photons from NR
	Int_t NumPhotonsER = NumPhotons - NumPhotonsNR; // Number of photons from ER
	cout << "There are " << NumPhotonsER << " photons from ER and " << NumPhotonsNR << " photons from NR" << endl;
	
	// Define functions for fast Sc fraction and for exp decay scintillation light
	TF1* FastProbFunc = new TF1 ("pdf for fast Sc","ROOT::Math::binomial_pdf(x,[0],[1])",0,NumPhotons+1);
	TF1 *ExpDecay     = new TF1 ("f_fast","exp(-x/[0])",0,50000*ns); // pdf of the decay
	
	// Simulate photons times for ER and NR
	cout << "ER=" << ER << " and NR=" << NR << endl;
	SimPhTimes(NumPhotonsER, ER, FastProbFunc, ExpDecay);
	SimPhTimes(NumPhotonsNR, NR, FastProbFunc, ExpDecay);
}
 
