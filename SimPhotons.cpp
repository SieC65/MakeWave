#include <iostream>

#include <TGraph.h>
#include "SystemOfUnits.h"

#include "SimPhotons.h"

using std::cout;
using std::endl;
using std::vector;
using CLHEP::ns;

SimPhotons::SimPhotons() {
	fMinPhotons = 100;
	fMaxPhotons = 10000;
	fTauFast    = 6*ns;
	fTauSlow    = 1500*ns;
	SetFastFrac (0.22, 0.75);
	gRandom->SetSeed(0);
	fFastProbFunc = new TF1 ("pdf for fast fraction","ROOT::Math::binomial_pdf(x,[0],[1])",0,100);
	fExpDecaySlow = new TF1 ("f_slow","exp(-x/[0])",0,100);
	fExpDecayFast = new TF1 ("f_fast","exp(-x/[0])",0,100);
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

vector <double> SimPhotons::SimulatePhotons(Int_t NumPhotons, Double_t FastFrac) {
	Int_t NumFast = round (NumPhotons * FastFrac);
	Int_t NumSlow = NumPhotons * (1 - FastFrac);
	SimulatePhotons (NumFast, NumSlow);
	return fSimPhotonTimes;
}

vector <double> SimPhotons::SimulatePhotons(Int_t NumFast, Int_t NumSlow) {

	// Set PDF of fast and slow scintillation decays
	Double_t XmaxSlow = 30*fTauSlow;
	fExpDecaySlow -> SetRange (0, XmaxSlow);
	fExpDecaySlow -> SetParameter (0, fTauSlow);

	Double_t XmaxFast = 30*fTauFast;
	fExpDecayFast -> SetRange (0, XmaxFast);
	fExpDecayFast -> SetParameter (0, fTauFast);

	fSimPhotonTimes.clear();
	// Get random tau for each photon
	for (Int_t phe = 0; phe < NumFast; phe++){
		fSimPhotonTimes.push_back(fExpDecayFast->GetRandom(0,XmaxFast));
	}
	for (Int_t phe = 0; phe < NumSlow; phe++) {
		fSimPhotonTimes.push_back(fExpDecaySlow->GetRandom(0,XmaxSlow));
	}

	return fSimPhotonTimes;
}

vector <double> SimPhotons::SimulatePhotons (Int_t NumPhotons, Option_t* type) {

	// Get fast fraction for this interaction type
	Double_t FastProb = 0;
	if (type == std::string("ER"))
		if (fFast_type == function)
			FastProb = fFastER_func -> Eval(NumPhotons);
		else
			FastProb = fFastER;
	else
	if (type == std::string("NR"))
		if (fFast_type == function)
			FastProb = fFastNR_func -> Eval(NumPhotons);
		else
			FastProb = fFastNR;
	else
		cout << "Please set \"NR\" or \"ER\" as interaction type" << endl;

	// Simulate number of fast photons
	fFastProbFunc -> SetRange (0, NumPhotons+1);
	fFastProbFunc -> SetParameter (0, FastProb);
	fFastProbFunc -> SetParameter (1, NumPhotons);
	Int_t NumFast = round (fFastProbFunc -> GetRandom (0, NumPhotons+1));


	// Simulate photons
	Int_t NumSlow = NumPhotons - NumFast;
	SimulatePhotons (NumFast, NumSlow);
	return fSimPhotonTimes;
}