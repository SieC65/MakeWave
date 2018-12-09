#include <iostream>

#include <TGraph.h>
#include "SystemOfUnits.h"

#include "SimPhotons.h"

using std::cout;
using std::endl;
using std::vector;
using CLHEP::ns;

SimPhotons::SimPhotons() {
	fMinPhotons = 140;
	fMaxPhotons = 2000;
	fTauFast    = 7*ns;
	fTauSlow    = 1600*ns;
	SetFastFrac (0.22, 0.75);
}

void SimPhotons::SetTau (Double_t TauFast, Double_t TauSlow) {
	fTauFast = TauFast;
	fTauFast = TauSlow;
}

void SimPhotons::SetFastFrac (TF1* FastER_func, TF1* FastNR_func) {
	fFast_type = function;
	fFastER_func = FastER_func;
	fFastNR_func = FastNR_func;
}

void SimPhotons::SetFastFrac (Double_t FastER, Double_t FastNR) {
	fFast_type = constant;
	fFastER = FastER;
	fFastNR = FastNR;	
}

void SimPhotons::SetDefFastFract () {
	fFast_type = function;
	fFastER_func = new TF1 ("Fast Sc Fraction for ER","[0]+[1]/x",fMinPhotons,fMaxPhotons);
	fFastER_func -> SetParameter (0, 0.178464);
    fFastER_func -> SetParameter (1, 46.705);
	fFastNR_func = new TF1 ("Fast Sc Fraction for NR","[0]+[1]/x",fMinPhotons,fMaxPhotons);
	fFastNR_func -> SetParameter (0, 0.723801);
    fFastNR_func -> SetParameter (1, -52.0528);
}

void SimPhotons::SimulatePhotons(Int_t NumPhotons, Double_t FastFrac) {
	Int_t NumFast = round (NumPhotons * FastFrac);
	Int_t NumSlow = NumPhotons * (1 - FastFrac);
	SimulatePhotons (NumFast, NumSlow);
}

void SimPhotons::SimulatePhotons(Int_t NumFast, Int_t NumSlow) {

	// Set PDF of fast and slow scintillation decays
	Double_t XmaxSlow = 3*fTauSlow;
	fExpDecaySlow = new TF1 ("f_slow","exp(-x/[0])",0,XmaxSlow);
	fExpDecaySlow -> SetParameter (0, fTauSlow);

	Double_t XmaxFast = 3*fTauFast;
	fExpDecayFast = new TF1 ("f_fast","exp(-x/[0])",0,XmaxFast);
	fExpDecayFast -> SetParameter (0, fTauFast);

	// Get random tau for each photon
	for (Int_t phe = 0; phe < NumFast; phe++){
		fSimPhotonTimes.push_back(fExpDecayFast->GetRandom(0,XmaxFast));
	}
	for (Int_t phe = 0; phe < NumSlow; phe++) {
		fSimPhotonTimes.push_back(fExpDecaySlow->GetRandom(0,XmaxSlow));
	}
}

void SimPhotons::SimulatePhotons (Int_t NumPhotons, Option_t* type) {
	if (type == std::string("ER"))
		fInterType = ER;
	else
	if (type == std::string("NR"))
		fInterType = NR;
	else
		cout << "Please set \"NR\" or \"ER\" as interaction type" << endl;
	Int_t NumFast = CalcNumFast (NumPhotons, fInterType);
	Int_t NumSlow = NumPhotons - NumFast;
	SimulatePhotons (NumFast, NumSlow);
}

Int_t SimPhotons::CalcNumFast(Int_t NumPhotons, InterType recoil) {

	// Calculate numbers of fast & slow photons for each group
	Double_t FastProb = 0;
	if      ((fFast_type == function) && (recoil == NR)){
		FastProb = fFastNR_func -> Eval(NumPhotons);
		cout << "NR case, function" << endl;
	}
	else if ((fFast_type == function) && (recoil == ER)){
		FastProb = fFastER_func -> Eval(NumPhotons);
		cout << "ER case, function" << endl;
	}
	else if ((fFast_type == constant) && (recoil == NR)){
		FastProb = fFastNR;
		cout << "NR case, const" << endl;
	}
	else if ((fFast_type == constant) && (recoil == ER)){
		FastProb = fFastER;
		cout << "ER case, const" << endl;
	}

	// Define function for fast Sc fraction
	TF1* FastProbFunc = new TF1 ("pdf for fast fraction","ROOT::Math::binomial_pdf(x,[0],[1])",0,NumPhotons+1);
	FastProbFunc -> SetRange (0, NumPhotons+1);
	FastProbFunc -> SetParameter (0, FastProb);
	FastProbFunc -> SetParameter (1, NumPhotons);

	// Simulate number of fast photons
	Int_t NumFast = round (FastProbFunc -> GetRandom (0, NumPhotons+1));
	Int_t NumSlow = NumPhotons - NumFast;

	return NumFast;
}