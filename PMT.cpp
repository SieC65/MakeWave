#include <iostream>
#include <TMath.h>
#include "PMT.h"

using namespace std;

//Simple constructor
PMT::PMT () {
	cout << "PMT object was created" << endl;
}

//Set PMT parameters
void PMT::SetPMTPar (Double_t QE, Double_t DPE, Double_t QE_1d, Double_t ArbAmpl_1d, Double_t GF_1d, Double_t AmplSigma, Double_t TOF_1d) {
	//Set values by user
	fQE			= QE;
	fDPE		= DPE;
	fQE_1d		= QE_1d;
	fArbAmpl_1d	= ArbAmpl_1d;
	fGF_1d		= GF_1d;
	fTOF_1d		= TOF_1d;
	//Calculate other values
	fQE_1d_ratio	= fArbAmpl_1d * fGF_1d * fQE_1d;	//auxiliary quantity
	fPphe_PC		= (fQE - fQE_1d_ratio)/(1 + fDPE - fQE_1d_ratio);	//Probability of photoeffect on PC
	fPphe_1d		= fGF_1d * fQE_1d;	//Probability of photoeffect on 1dyn for photons passed PC
	fAmplSigma_1d	= AmplSigma * fArbAmpl_1d;	//Sigma for amplitude from 1dyn
}

Double_t PMT::GetPphe_PC () {
	return fPphe_PC;
}

Double_t PMT::GetPphe_1d () {
	return fPphe_1d;
}

Double_t PMT::GetArbAmpl_1d () {
	return fArbAmpl_1d;
}

Double_t PMT::GetAmplSigma_1d () {
	return fAmplSigma_1d;
}

Double_t PMT::GetDPE() {
	return fDPE;
}

Double_t PMT::GetTOF_1d() {
	return fTOF_1d;
}