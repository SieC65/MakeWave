#include <iostream>

#include <TGraph.h>
#include "SystemOfUnits.h"

#include "MakeWave.h"

using std::cout;
using std::endl;
using std::vector;
using CLHEP::ns;

//using CLHEP::mV;

// Simple constructor
MakeWave::MakeWave () {
	PhotonPulse = 0;
	DarkPulse   = 0;
	fPMT        = 0;
	//cout << "MakeWave object was created" << endl;
}

// Set PMT
void MakeWave::SetPMT (RED::PMT* PMT) {
	fPMT = PMT;
	//cout << "PMT was set" << endl;
}

// Set OutWave parameters
void MakeWave::SetOutWave (Double_t Period, Double_t Gain, Int_t NumSamples, Double_t Delay) {
	fPeriod     = Period;
	fGain       = Gain;
	fNumSamples = NumSamples;
	fDelay      = Delay;
	//cout << "OutWave parameters were set" << endl;
}

// Set sequence of photon times
void MakeWave::SetPhotonTimes (vector <double>* PhotonTimes) {
	fPhotonTimes = PhotonTimes;
	//cout << "Sequence of SPE's arrival times was set" << endl;
}

// Creating OutWave
void MakeWave::CreateOutWave () {
	//cout << "Creating OutWave..." << endl;

	fOutWave.clear();              //  Clear vector OutWave
	fOutWave.resize (fNumSamples, 0); // Resize vector OutWave
	Char_t NumPhe = 0;
	
	// Create PulseArray vector
	if (!PhotonPulse) {
		PhotonPulse = new RED::PMT::PulseArray;
	}
	else PhotonPulse->clear();

	// Create DarkPulseArray vector
	if (!DarkPulse) {
		DarkPulse = new RED::PMT::PulseArray;
	}
	else DarkPulse->clear();

	/// ADD SPE FROM PHOTONS
	//cout << "adding photons" << endl;
	// Generate vector Pulse of PulseArray type
	for (unsigned int i = 0; i < fPhotonTimes->size(); i++) {
		NumPhe = fPMT->OnePhoton (&(fPhotonTimes->at(i)), *PhotonPulse, true);
	}
	AddPulseArray (PhotonPulse);

	/// ADD DARK COUNTS
	//cout << "adding dark" << endl;
	// Generate vector DarkPulse of PulseArray type
	fPMT->GenDCR (fDelay - (fPMT->GetXmax() - fPMT->GetXmin()), fDelay + fNumSamples * fPeriod, *DarkPulse);
	AddPulseArray (DarkPulse);
	//cout << "OutWave was created" << endl;
}

// Add PulseArray vector to OutWave
void MakeWave::AddPulseArray (RED::PMT::PulseArray *Pulses) {
	Double_t SampleTime   = 0; // Time of sample from "0" of OutWave
	Double_t PulseTime    = 0; // Time from "0" of OutWave to "0" of SPE shape
	Double_t PulseAmpl    = 0; // Amplitude of SPE shape
	Int_t StartSample     = 0; // First sample in SPE domain
	Int_t FinishSample    = 0; // Last sample in SPE domain
	
	// Go along all pulses in PulseArray vector DarkPulse and add them to OutWave
	for (unsigned int i = 0; i < Pulses->size(); i++) {
		// Calculate left and right samples including SPE
		StartSample     =  ceil( (((*Pulses)[i]).fTime - fDelay + fPMT->GetXmin()) / fPeriod );
		FinishSample    = floor( (((*Pulses)[i]).fTime - fDelay + fPMT->GetXmax()) / fPeriod );
		// Get delay time from "0" of OutWave to "0" of SPE shape
		PulseTime       = ((*Pulses)[i]).fTime - fDelay;
		// Get amplitude of SPE shape
		PulseAmpl       = ((*Pulses)[i]).fAmpl;
		// Limit edges
		if (StartSample < 0)
			StartSample = 0;
		if (FinishSample > fNumSamples - 1)
			FinishSample = fNumSamples - 1;
		// Add SPE to OutWave
		for (int s = StartSample; s <= FinishSample; s++) {
			SampleTime = s*fPeriod;
			fOutWave[s] += PulseAmpl/fGain * fPMT->Eval(SampleTime - PulseTime);
		}
	}
}

// Draws histograms
void MakeWave::test(RED::PMT_R11410 &) {
	// Histograms
	TH1F* TimeHist;      // Delay time from photon hit to pulse
	TH1F* AmplHist;      // Amplitude of pulse from photon	
	TH1F* DarkTimeHist;  // Time of dark pulses
	TH1F* DarkAmplHist;  // Amplitude of dark pulses
	TH1F *NumPheHist;    // Number of photoelectron created by 1 photon:
	                     // -1 or -2: 1 or 2 phe in 1dyn
	                     // +1 or +2: 1 or 2 phe in PC
	TH1F *PulseAreaHist; // Area under pulse SPE (DPE)
}

// Print OutWave
void MakeWave::PrintOutWave() {
	cout << "Printing OutWave..." << endl;
	for (int i = 0; i < fNumSamples; i++) {
		cout << fOutWave[i] << "\t(t = " << (fDelay + i*fPeriod)/ns << " ns)" << endl;
	}
	cout << "OutWave was printed" << endl;
}

// Draw OutWave
void MakeWave::DrawOutWave () {
	cout << "Drawing OutWave..." << endl;
	TGraph *g1 = new TGraph(fNumSamples);
	for (int i = 0; i < fNumSamples; i++) {
		g1->SetPoint(i, (fDelay + i*fPeriod)/ns, fOutWave[i]);
	}
	g1->SetTitle("Waveform");
	g1->Draw();
	cout << "OutWave was drawn" << endl;
}

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
