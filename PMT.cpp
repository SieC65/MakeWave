#include <iostream>
#include <TMath.h>
#include "PMT.h"

using namespace std;

//Simple constructor
PMT::PMT () {
	cout << "PMT object was created" << endl;
}

//Set PMT parameters
void PMT::SetPMTPar (const PMTPar &PMTX, Double_t AmplSigma) {
	//Set values by user
	fPMTX = PMTX;
	//Calculate other values
	fQE_1d_ratio	= fPMTX.ArbAmpl_1d * fPMTX.GF_1d * fPMTX.QE_1d;	//auxiliary quantity
	fPphe_PC		= (fPMTX.QE - fQE_1d_ratio)/(1 + fPMTX.DPE - fQE_1d_ratio);	//Probability of photoeffect on PC
	fPphe_1d		= fPMTX.GF_1d * fPMTX.QE_1d;	//Probability of photoeffect on 1dyn for photons passed PC
	fAmplSigma_1d	= AmplSigma * fPMTX.ArbAmpl_1d;	//Sigma for amplitude from 1dyn
}