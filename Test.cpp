#include <iostream>

#include <TGraph.h>
#include "SystemOfUnits.h"

#include "Test.h"

using std::cout;
using std::endl;
using std::vector;
using CLHEP::ns;

// MAKEPHOTONS

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
	if      ((fFast_type == function) && (recoil == ER)){
		FastProb = fFastER_func -> Eval(NumPhotons);
		cout << "ER case, function" << endl;
	}
	else if ((fFast_type == function) && (recoil == NR)){
		FastProb = fFastNR_func -> Eval(NumPhotons);
		cout << "NR case, function" << endl;
	}
	else if ((fFast_type == constant) && (recoil == ER)){
		FastProb = fFastER;
		cout << "ER case, const" << endl;
	}
	else if ((fFast_type == constant) && (recoil == NR)){
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
	for (Int_t phe = 0; phe < NumFast; phe++){
		fSimPhotonTimes.push_back(ExpDecay->GetRandom(0,500*ns));
	}
	ExpDecay -> SetRange (0, 50000*ns);
	ExpDecay -> SetParameter (0, fTauSlow);
	for (Int_t phe = 0; phe < NumSlow; phe++) {
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
 
// MAKETEST

MakeTest::MakeTest() {
	fTimesVec = new vector <Double_t>;
	fMeanOutWave = new vector <Double_t>;
	cout << "MakeTest object created" << endl;
}

void MakeTest::SetPhotonsTimes (Double_t *TimesArr) {
	Int_t TimesArrSize = sizeof(TimesArr)/sizeof(TimesArr[0]);
	(*fTimesVec).resize (TimesArrSize);
	for (int i = 0; i < TimesArrSize; i++) {
		(*fTimesVec)[i] = TimesArr[i]*ns;
	}
}

void MakeTest::AverageOW (Int_t DebugN, MakeWave* MakeWaveObj) {

	//Double_t Period     = MakeWaveObj -> GetPeriod();
	//Double_t Gain       = MakeWaveObj -> GetGain();
	Double_t NumSamples = MakeWaveObj -> GetNumSamples();
	//Double_t Delay      = MakeWaveObj -> GetDelay();

	cout << "Create vector for mean waveform" << endl;
	fMeanOutWave->resize (NumSamples);
	Bool_t ShowAllNum = false; // Report status of each OutWave performing
	                           // If false - only 10 steps will be showed
	// Create OutWave
	if (DebugN < 100) {
		ShowAllNum = true;
	}
	cout << "Get into cycle" << endl;
	for (int i = 1; i <= DebugN; i++){
		if (ShowAllNum)
			cout << "Creating " << i << " waveform" << endl;
		else{
			if (!(i % int(floor(DebugN/10)))){
				cout << "Creating " << i << " waveform" << endl;
			}
		}
		MakeWaveObj->CreateOutWave();
		for (int k = 0; k < NumSamples; k++){
			(*fMeanOutWave)[k] += MakeWaveObj->GetOutWave()[k];
		}
	}
}

void MakeTest::DrawOW (MakeWave* MakeWaveObj) {
	Double_t Period     = MakeWaveObj -> GetPeriod();
	//Double_t Gain       = MakeWaveObj -> GetGain();
	Double_t NumSamples = MakeWaveObj -> GetNumSamples();
	Double_t Delay      = MakeWaveObj -> GetDelay();

	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->SetTitle("OutWave");
	TGraph *g2  = new TGraph (NumSamples);
	for (int i = 0; i < NumSamples; i++){
		g2->SetPoint (i, (Delay + i*Period)/ns, (*fMeanOutWave)[i]);
		cout << "point num " << i << "at time= " << (Delay + i*Period)/ns << " ns ";
		cout << "with ampl= " << (*fMeanOutWave)[i] << " ADC-units was added" <<endl;
	}
	g2->SetTitle("OutWave;Time, ns;Amplitude, ADC units");
	g2->Draw();
}

void MakeTest::RandPhotonTimes (Int_t number, Double_t leftEdge, Double_t rightEdge) {
	fTimesVec -> resize (number);
	TRandom3 rand;
	rand.SetSeed(0);
	for (int i=0; i<number; i++){
		(*fTimesVec)[i] = leftEdge + rand.Rndm() * (rightEdge - leftEdge);
	}
	sort((*fTimesVec).begin(), (*fTimesVec).end());
}

void MakeTest::PrintTimesVec () {
	cout << "Arrival times:" << endl;
	for (unsigned int i = 0; i < fTimesVec->size(); i++){
		cout << ((*fTimesVec)[i])/ns << " ns" << endl;
	}
}