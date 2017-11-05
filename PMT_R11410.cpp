#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TMath.h>

#include "PMT_R11410.hh"
#include "SystemOfUnits.h"

using std::cout;
using std::endl;
namespace CLHEP {
	static const double mV = 1.e-3*volt;
}
using CLHEP::ns;
using CLHEP::mV;

namespace RED
{	
	PMT_R11410::PMT_R11410() {
		cout << "PMT_R11410 object was created" << endl;
	}
	
	//	Set PMT Hamamatsu R11410-20 parameters.
	//	Some values are taken from articles C.H.Faham et.al. "Measurements of wavelength-dependent
	//	double photoelectron emission..." //JINST 2015 and D.Yu.Akimov et.al. "Performance of
	//	Hamamatsu R11410-20 PMTs..." //JINST 2016
	void PMT_R11410::SetDefaults() {		
		// Set independent parameters of PMT
		fQE			= 0.3;		// QE (standard, i.e. with many phe)
		fDPE		= 0.225;	// Double Photoelectron Emission P(2phe)/(P(2phe)+P(1phe))
		fQE_1d		= 0.105;	// QE for 1dyn (only 1phe), prob.
		fGain_PC_1d	= 13;		// Amplification on first gap (PC-1dyn)
		fGF_1d		= 0.1;		// Average geom. prob. for a rndm photon from PC to hit 1st dyn
		fArea_mean	= 10*mV*ns;	// SPE pulse area
		fArea_sigma	= 1*mV*ns;	//
		fTOFe_PC_1d	= 6*ns;		// ToF e- from PC to 1dyn
		fTOFe_mean	= 30*ns;	// ToF e- from PC to anode
		fTOFe_sigma	= 3*ns;		//
		fDCR		= 100;		// Dark count rate, Hz
		fAP_cont	= 0;		//
		fAP_peak	= 0;		// Afterpulsing probability
		
		// Set shape as gaussian
		fMode		= kModeF1;
		fShape.f1	= new TF1("SPE","gaus(0)",-5*ns,5*ns);
		fShape.f1->SetParameter(0, 1);		//Amplitude of gaussian ( = 1 )
		fShape.f1->SetParameter(1, 0);		//Center
		fShape.f1->SetParameter(2, 1*ns);	//Sigma
		
		cout << endl << "Default PMT parameters (QE, DPE, SPE Area, TOF(e), DCR, AP, GF(1d), Gain(PC->1d)) were set" << endl;
		
		// Calculate all dependent parameters
		CalculateParams();
	}

	void PMT_R11410::SetParameters(	double_t QE, double_t DPE, double_t QE_1d, double_t Gain_PC_1d, 
									double_t GF_1d, double_t Area_mean, double_t Area_sigma, 
									double_t TOFe_PC_1d, double_t TOFe_mean, double_t TOFe_sigma, 
									double_t DCR, double_t AP_cont, double_t AP_peak) {
		fQE			= QE;
		fDPE		= DPE;
		fQE_1d		= QE_1d;
		fGain_PC_1d	= Gain_PC_1d;
		fGF_1d		= GF_1d;
		fArea_mean	= Area_mean;
		fArea_sigma	= Area_sigma;
		fTOFe_PC_1d	= TOFe_PC_1d;
		fTOFe_mean	= TOFe_mean;
		fTOFe_sigma	= TOFe_sigma;
		fDCR		= DCR;
		fAP_cont	= AP_cont;
		fAP_peak	= AP_peak;
		cout << endl << "User-defined PMT parameters (QE, DPE, SPE Area, TOF(e), DCR, AP, GF(1d), Gain(PC->1d)) were set" << endl;
		CalculateParams();
	}	
	
	void PMT_R11410::SetShape(TF1* Shape) {
		fMode = kModeF1;
		fShape.f1	= Shape;
		cout << endl << "User-defined PMT Shape was set as TF1" << endl;
	}
	
	void PMT_R11410::SetShape(TSpline* Shape) {
		fMode = kModeSpline;
		fShape.f2	= Shape;
		cout << endl << "User-defined PMT Shape was set as TSpline" << endl;
	}
	
	void PMT_R11410::DrawShape(const char* title) {
		TCanvas *c5 = new TCanvas();
		c5->SetTitle("SPE Shape");
		c5->cd();
		if (fMode==1) {
			fShape.f1->SetTitle(title);
			fShape.f1->Draw();
		} else if (fMode==2) {
			fShape.f2->SetTitle(title);
			fShape.f2->Draw();
		}
		//c5->WaitPrimitive();
	}
	
	void PMT_R11410::CalculateParams() {
		// Get probabilities of interaction
		fProb_C		= fQE/(1+fDPE) - fQE_1d;	// Prob. of int. in PC
		fProb_C1 	= (1-fDPE) * fProb_C;		// Prob. of 1phe in PC
		fProb_C2 	= fDPE * fProb_C;			// Prob. of 2phe in PC
		fProb_1d 	= fQE_1d * fGF_1d * (1-fProb_C);	// Full prob. of int. in 1d
		
		// Get Area of SPE from 1d (assume that the width is the same as one from PC)
		fArea_1d_mean 	= fArea_mean / fGain_PC_1d;
		fArea_1d_sigma 	= fArea_sigma / fGain_PC_1d;
		
		// Get Amplitude from Area
		if (fMode==1)
			fAmplitude_mean = fArea_mean / fShape.f1->Integral(GetXmin (), GetXmax ()) ;
		else if (fMode==2)
			fAmplitude_mean = 1 ;
		else
			fAmplitude_mean = 0;
		fAmplitude_sigma 	= fAmplitude_mean * (fArea_sigma/fArea_mean) ;
		
		// Get TOF(e-) from 1d to Anode
		fTOFe_1d_mean 	= fTOFe_mean - fTOFe_PC_1d;
		fTOFe_1d_sigma 	= 0.5 * fTOFe_sigma;		// May be wrong
		
		cout << "Additional PMT parameters (interaction probabilities, Amplitude od SPE, Area and TOF(e) from 1dyn) were calculated" << endl;
	}

	Double_t PMT_R11410::GetXmax () const {
		if (fMode==1)
			return(fShape.f1->GetXmax());
		else if (fMode==2)
			return(fShape.f2->GetXmax());
		else
			return 0;
	}

	Double_t PMT_R11410::GetXmin () const {
		if (fMode==1)
			return(fShape.f1->GetXmin());
		else if (fMode==2)
			return(fShape.f2->GetXmin());
		else
			return 0;
	}

	Double_t PMT_R11410::GetMaximum () const {
		if (fMode==1)
			return(fShape.f1->GetMaximum());
		else if (fMode==2) {
			cout << "The GetMaximum function is not yet implemented for spline" << endl;
			return 0;
		}
		else
			return 0;
	}

	Double_t PMT_R11410::GetMinimum () const {
		if (fMode==1)
			return(fShape.f1->GetMinimum());
		else if (fMode==2) {
			cout << "The GetMinimum function is not yet implemented for spline" << endl;
			return 0;
		}
		else
			return 0;
	}

	Double_t PMT_R11410::Eval (Double_t t) const {
		if (fMode==1)
			return(fShape.f1->Eval(t));
		else if (fMode==2)
			return(fShape.f2->Eval(t));
		else
			return 0;
	}
	
	int PMT_R11410::Begin(PulseArray &electrons) {
		cout << "Function 'Begin' is not yet implemented" << endl;
		return(0);
	}

	int PMT_R11410::End(PulseArray &electrons) {
		cout << "Function 'End' is not yet implemented" << endl;
		return(0);
	}

	// Convert photon to pulses. time - time of arr. photon
	int PMT_R11410::OnePhoton(double_t time, PulseArray &electrons, bool fDebug) {
		Pulse OnePulse;		// Temporary variable for saving each SPE
		Double_t RND = 0; 	// Random in range 0..1 defining interact. type
			// We can divide range 0..1 into 4 bands:
			// 0 .. fProb_C1 												- it was SPE from PC
			// fProb_C1 .. fProb_C1 + fProb_C2 								- it was DPE from PC
			// (fProb_C1 + fProb_C2) .. (1fProb_C1 + fProb_C2) + fProb_1d 	- it was SPE from 1dyn
			// (1fProb_C1 + fProb_C2 + fProb_1d) .. 1 						- it was no interaction
		fRND.SetSeed(0);
		RND = fRND.Rndm();
		Int_t InteractType = 0;	//0-no int, 1 - 1phe in PC, 2 - 2phe in PC, -1 - 1phe in 1d
		
if (fDebug) cout << "At " << time/ns << " ns ";

		// If photon created 1 phe in PC
		if (RND < fProb_C1) {
if (fDebug) cout << "photon created 1 phe in PC" << endl;
			InteractType = 1;
			OnePulse.fAmplitude	= fRND.Gaus (fAmplitude_mean, fAmplitude_sigma);
			OnePulse.fTime		= fRND.Gaus (fTOFe_mean, fTOFe_sigma) + time;
if (fDebug) cout << "with amplitude = " << OnePulse.fAmplitude/mV << " mV and time from '0' of OutWave =" << OnePulse.fTime/ns << " ns" << endl;
			electrons.push_back(OnePulse);		
		}
		
		// If photon created 2 phe in PC
		else if (RND < fProb_C1 + fProb_C2) {
if (fDebug) cout << "photon created 2 phe in PC" << endl;
			InteractType = 2;
			for (int i = 0; i < 2; i++) {
				OnePulse.fAmplitude	= fRND.Gaus (fAmplitude_mean, fAmplitude_sigma);
				OnePulse.fTime		= fRND.Gaus (fTOFe_mean, fTOFe_sigma) + time;
if (fDebug) cout << "with amplitude = " << OnePulse.fAmplitude/mV << " mV and time from '0' of OutWave =" << OnePulse.fTime/ns << " ns" << endl;
				electrons.push_back(OnePulse);
			}			
		}
			
		// If photon created phe in 1dyn
		else if (RND < fProb_C1 + fProb_C2 + fProb_1d) {
if (fDebug) cout << "photon created phe in 1dyn" << endl;
			InteractType = -1;
			OnePulse.fAmplitude	= fRND.Gaus (fAmplitude_mean/fGain_PC_1d, fAmplitude_sigma);
			OnePulse.fTime		= fRND.Gaus (fTOFe_1d_mean, fTOFe_1d_sigma) + time;
if (fDebug) cout << "with amplitude = " << OnePulse.fAmplitude/mV << " mV and time from '0' of OutWave =" << OnePulse.fTime/ns << " ns" << endl;
			electrons.push_back(OnePulse);
		}
		
		// Else photon didn't interact with PMT, nothing to do
		else {
if (fDebug) cout << "photon didn't interact with PMT" << endl;
			InteractType = 0;
		}
		
		return(InteractType);
	}
	
	void PMT_R11410::Clear(Option_t *option) { 
		cout << "Function 'Clear' is not yet implemented" << endl;
	}

	void PMT_R11410::Print(Option_t *option) const {
		cout << endl << "Printing PMT parameters..." << endl;
		cout << "Defined parameters:" << endl;
		cout << "\t SPE shape type \t=\t";
		switch (fMode) {
			case 0:
				cout << "None (still can cause fault)" << endl;
				break;
			case 1:
				cout << "TF1" << endl;
				break;
			case 2:
				cout << "TSpline (still can cause fault)" << endl;
				break;
		}
		cout << "\t SPE Xmin\t\t=\t" << GetXmin()/ns << " ns \t\t//SPE shape Xmin" << endl;
		cout << "\t SPE Xmax\t\t=\t" << GetXmax()/ns << " ns \t\t//SPE shape Xmax" << endl;
		cout << "\t SPE Area mean\t\t=\t" << fArea_mean/(mV*ns) << " mV*ns \t//SPE Area, mean" << endl;
		cout << "\t SPE Area sigma\t\t=\t" << fArea_sigma/(mV*ns) << " mV*ns \t//SPE Area, sigma" << endl;
		cout << "\t QE of PC\t\t=\t" << fQE*100 << " % \t\t//Quantum Efficiency of Photocathode" << endl;
		cout << "\t QE of 1d\t\t=\t" << fQE_1d*100 << " % \t\t//Quantum Efficiency of 1st dynode" << endl;
		cout << "\t DPE\t\t\t=\t" << fDPE*100 << " % \t\t//Double Photoelectron Emission probability" << endl;
		cout << "\t GF of 1d\t\t=\t" << fGF_1d*100 << " % \t\t//Geometric factor (or probability for a random photon from PC to hit 1st dynode)" << endl;
		cout << "\t PC->1d gain\t\t=\t" << fGain_PC_1d << "\t\t//Amplification in the Photocathode - 1st dynode gap" << endl;
		cout << "\t DCR\t\t\t=\t" << fDCR << " Hz \t\t//Dark count rate" << endl;
		cout << "\t TOF PC->An mean\t=\t" << fTOFe_mean/ns << " ns \t\t//Time of Flight from Photocathode to anode, mean" << endl;
		cout << "\t TOF PC->An sigma\t=\t" << fTOFe_sigma/ns << " ns \t\t//Time of Flight from Photocathode to anode, sigma" << endl;
		cout << "\t TOF PC->1d mean\t=\t" << fTOFe_PC_1d/ns << " ns \t\t//Time of Flight from Photocathode to 1st dynode, mean" << endl;
		cout << "\t AP cont\t\t=\t" << fAP_cont << "\t\t//Afterpulses ?" << endl;
		cout << "\t AP peak\t\t=\t" << fAP_peak << " \t\t//Afterpulses ?" << endl;
		
		cout << "Calculated parameters:" << endl;
		cout << std::setiosflags(std::ios::fixed) << std::setprecision(2);
		cout << "\t Total PC probability \t=\t" << fProb_C*100 << " % \t//Total interaction probability on Photocathode" << endl;
		cout << "\t 1phe PC probability \t=\t" << fProb_C1*100 << " % \t//1phe interaction probability on Photocathode" << endl;
		cout << "\t 2phe PC probability \t=\t" << fProb_C2*100 << " % \t\t//2phe interaction probability on Photocathode" << endl;
		cout << "\t Total 1d probability \t=\t" << fProb_1d*100 << " % \t\t//Total interaction probability on 1st dynode" << endl;
		cout << std::resetiosflags(std::ios::fixed) << std::setprecision(6);
		cout << "\t SPE 1d Area mean \t=\t" << fArea_1d_mean/(mV*ns) << " mV*ns \t//1st dynode SPE Area, mean" << endl;
		cout << "\t SPE 1d Area sigma \t=\t" << fArea_1d_sigma/(mV*ns) <<" mV*ns \t//1st dynode SPE Area, sigma" << endl;
		cout << "\t TOF 1d->An mean \t=\t" << fTOFe_1d_mean/ns << " ns \t\t//Time of Flight from 1st dynode to anode, mean" << endl;
		cout << "\t TOF 1d->An sigma \t=\t" << fTOFe_1d_sigma/ns << " ns \t\t//Time of Flight from 1st dynode to anode, sigma" << endl;
		if (fMode == 1) {
			cout << "\t SPE Shape Area \t=\t" << fShape.f1->Integral(GetXmin (), GetXmax ())/(ns) << " arb.un.*ns \t//SPE shape Area, integrated from Xmin to Xmax by ROOT standard method" << endl;
			cout << "\t SPE Shape MinValue\t=\t" << GetMinimum() << " arb.un. \t//Minimum value of SPE shape in range Xmin .. Xmax" << endl;
			cout << "\t SPE Shape MaxValue\t=\t" << GetMaximum() << " arb.un. \t//Maximum value of SPE shape in range Xmin .. Xmax" << endl;
		}
	}
}