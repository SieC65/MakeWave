#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TMath.h>

#include "PMT_R11410.hh"

using std::cout;
using std::endl;

namespace RED
{
	PMT_R11410::PMT_R11410() {
		cout << "PMT_R11410 object was created" << endl;
	}

	// Set PMT Hamamatsu R11410-20 parameters.
	// Some values are taken from articles C.H.Faham et.al. "Measurements of wavelength-dependent
	// double photoelectron emission..." //JINST 2015 and D.Yu.Akimov et.al. "Performance of
	// Hamamatsu R11410-20 PMTs..." //JINST 2016
	void PMT_R11410::SetDefaults() {        
		// Set independent parameters of PMT
		fQE         = 0.3;      // QE (standard, i.e. with many phe)
		fDPE_PC     = 0.225;    // Double Photoelectron Emission probability for PC
		fQE_1d      = 0.105;    // QE for 1dyn
		fDPE_1d     = 0;        // Double Photoelectron Emission probability for 1dyn
		fGain_PC_1d = 13;       // Amplification on first gap (PC-1dyn)
		fGF_1d      = 0.1;      // Average geom. prob. for a rndm photon from PC to hit 1dyn
		fArea_mean  = 10*mV*ns; // SPE pulse area
		fArea_sigma = 1*mV*ns;  // SPE pulse area error
		fTOFe_PC_1d = 6*ns;     // ToF e- from PC to 1dyn
		fTOFe_mean  = 30*ns;    // ToF e- from PC to anode
		fTOFe_sigma = 3*ns;     // ToF e- from PC to anode error
		fDCR        = 100;      // Dark count rate, Hz
		fAP_cont    = 0;        // ???
		fAP_peak    = 0;        // Afterpulsing probability
		
		// Set shape as gaussian
		fMode     = kModeF1;
		fShape.func = new TF1("SPE","gaus(0)",-5*ns,5*ns);
		fShape.func->SetParameter(0, 1);    // Amplitude of gaussian (1)
		fShape.func->SetParameter(1, 0*ns); // Center of gaussian (0 ns)
		fShape.func->SetParameter(2, 1*ns); // Sigma of gaussian (1 ns)
		
		cout << endl;
		cout << "Default PMT parameters ";
		cout << "(QE, DPE, SPE Area, TOF(e), DCR, AP, GF(1d), Gain(PC->1d)) ";
		cout << "were set" << endl;
		
		// Calculate all dependent parameters
		CalculateParams();
	}

	void PMT_R11410::SetParams (double QE, double Area_mean, double DCR, double AP_cont) {
		fQE         = QE;
		fArea_mean  = Area_mean;
		fDCR        = DCR;
		fAP_cont    = AP_cont;
		cout << endl << "User-defined PMT parameters ";
		cout << "QE, SPE Area, DCR, AP were set as follows:" << endl;
		cout << "QE       = " << fQE << endl;
		cout << "SPE Area = " << fArea_mean/(mV*ns) << " mV*ns" << endl;
		cout << "DCR      = " << fDCR << " Hz" << endl;
		cout << "AP       = " << fAP_cont << endl;
		CalculateParams();
	}

	void PMT_R11410::SetDPE_PC (Double_t DPE_PC) {
		fDPE_PC = DPE_PC;
		cout << "DPE_PC was set as " << fDPE_PC << endl;
	}

	void PMT_R11410::SetDPE_1d (Double_t DPE_1d) {
	fDPE_1d = DPE_1d;
	cout << "DPE_1d was set as " << fDPE_1d << endl;
}

	void PMT_R11410::SetQE_1d (Double_t QE_1d) {
		fQE_1d = QE_1d;
		cout << "QE at 1st dyn was set as " << fQE_1d << endl;
	}

	void PMT_R11410::SetGain_PC_1d (Double_t Gain_PC_1d) {
		fGain_PC_1d = Gain_PC_1d;
		cout << "Gain(PC->1d) was set as " << fGain_PC_1d << endl;
	}

	void PMT_R11410::SetGF_1d (Double_t GF_1d) {
		fGF_1d = GF_1d;
		cout << "Geom factor of 1st dyn was set as " << fGF_1d << endl;
	}

	void PMT_R11410::SetArea_sigma (Double_t Area_sigma) {
		fArea_sigma = Area_sigma;
		cout << "SPE area sigma was set as " << fArea_sigma/(mV*ns) << " mV*ns" << endl;
	}

	void PMT_R11410::SetTOFe_PC_1d (Double_t TOFe_PC_1d) {
		fTOFe_PC_1d = TOFe_PC_1d;
		cout << "TOF(e) PC->1d was set as " << fTOFe_PC_1d << endl;
	}

	void PMT_R11410::SetTOFe_mean (Double_t TOFe_mean) {
		fTOFe_mean = TOFe_mean;
		cout << "TOF(e) was set as " << fTOFe_mean << endl;
	}

	void PMT_R11410::SetTOFe_sigma (Double_t TOFe_sigma) {
		fTOFe_sigma = TOFe_sigma;
		cout << "TOF(e) sigma was set as " << fTOFe_sigma << endl;
	}

	void PMT_R11410::SetAP_peak (Double_t AP_peak) {
		fAP_peak = AP_peak;
		cout << "AP peak was set as " << fAP_peak << endl;
	}

	void PMT_R11410::SetShape (TF1* Shape) {
		fMode = kModeF1;
		fShape.func   = Shape;
		cout << endl << "User-defined PMT Shape was set as TF1" << endl;
	}

	void PMT_R11410::SetShape (TSpline* Shape) {
		fMode = kModeSpline;
		fShape.spline   = Shape;
		cout << endl << "User-defined PMT Shape was set as TSpline" << endl;
	}

	void PMT_R11410::DrawShape (const char* title) {
		TCanvas *c5 = new TCanvas();
		c5->SetTitle("SPE Shape");
		c5->cd();
		switch (fMode) {
			case kModeNone:
				break;
			case kModeF1 :
				fShape.func->SetTitle (title);
				fShape.func->Draw();
				break;
			case kModeSpline :
				fShape.spline->SetTitle (title);
				fShape.spline->Draw();
				break;
		}
		//c5->WaitPrimitive();
	}

	void PMT_R11410::CalculateParams() {
		// Get probabilities of interaction
		double_t fQE_1d_ratio = 1./fGain_PC_1d * fGF_1d * fQE_1d; // Number of phe (given to PC, as if these phe 
		                                            // were born in PC) per 1 photon passed through PC:
		                                            // fQE = fProb_C * (1 + fDPE_PC) + (1 - fProb_C) * fQE_1d_ratio
		fProb_C   = (fQE - fQE_1d_ratio) / (1 + fDPE_PC - fQE_1d_ratio); //  Prob. of int. in PC
		fProb_1d  = (1 - fProb_C) * fGF_1d * fQE_1d / (1 + fDPE_1d); // Full prob. of int. in 1d
		fProb_C1  = (1 - fDPE_PC) * fProb_C;  // Prob. of 1phe in PC
		fProb_C2  = fDPE_PC       * fProb_C;  // Prob. of 2phe in PC
		fProb_1d1 = (1 - fDPE_1d) * fProb_1d; // Full prob. of 1phe in 1d
		fProb_1d2 = fDPE_1d       * fProb_1d; // Full prob. of 2phe in 1d

		// Get Area of SPE from 1d (assume that the width is the same as one from PC
		// and the area is proportional to the amplitude)
		fArea_1d_mean   = fArea_mean  / fGain_PC_1d;
		fArea_1d_sigma  = fArea_sigma / fGain_PC_1d;

		// Get Amplitude from Area
		switch (fMode) {			case kModeF1 :   
				fAmpl_mean = fArea_mean / fShape.func->Integral(GetXmin (), GetXmax ()) ;
				break;
			case kModeSpline :
				fAmpl_mean = fArea_mean / SplineIntegral (fShape.spline, GetXmin (), GetXmax (), 1000) ;
				break;
			default:
				fAmpl_mean = 0;
				break;
		}
		fAmpl_sigma    = fAmpl_mean * (fArea_sigma/fArea_mean) ;

		// Get TOF(e-) from 1d to Anode
		fTOFe_1d_mean   = fTOFe_mean - fTOFe_PC_1d;
		fTOFe_1d_sigma  = 0.5 * fTOFe_sigma; // May be wrong
		
		cout << "Additional PMT parameters ";
		cout << "(interaction probabilities, Amplitude od SPE, Area and TOF(e) from 1dyn) ";
		cout << "were calculated" << endl;
	}

	Double_t PMT_R11410::GetXmax () const {
		switch (fMode) {			case kModeF1 :   
				return(fShape.func->GetXmax());
				break;
			case kModeSpline :
				return(fShape.spline->GetXmax());
				break;
			default:
				return 0;
				break;
		}
	}

	Double_t PMT_R11410::GetXmin () const {
		switch (fMode) {			case kModeF1 :   
				return(fShape.func->GetXmin());
				break;
			case kModeSpline :
				return(fShape.spline->GetXmin());
				break;
			default:
				return 0;
				break;
		}
	}

	Double_t PMT_R11410::GetYmax () const {
		switch (fMode) {			case kModeF1 :   
				return(fShape.func->GetMaximum());
				break;
			case kModeSpline :
				cout << "The GetYmax function is not yet implemented for spline" << endl;
				return 0;
				break;
			default:
				return 0;
				break;
			}
	}

	Double_t PMT_R11410::GetYmin () const {
		switch (fMode) {			case kModeF1 :   
				return (fShape.func->GetMinimum());
				break;
			case kModeSpline :
				cout << "The GetYmin function is not yet implemented for spline" << endl;
				return 0;
				break;
			default:
				return 0;
				break;
			}
	}

	Double_t PMT_R11410::Eval (Double_t t) const {
		switch (fMode) {			case kModeF1 :   
				return (fShape.func->Eval(t));
				break;
			case kModeSpline :
				return (fShape.spline->Eval(t));
				break;
			default:
				return 0;
				break;
		}
	}

	Double_t PMT_R11410::SplineIntegral (TSpline* spline, Double_t xmin, Double_t xmax, Int_t nbins) {
		Double_t Integral = 0;
		Double_t BinWidth = (xmax - xmin) / nbins;
		for (int i=0; i<nbins; i++) {
			Integral += BinWidth * spline->Eval(xmin + (i+0.5)*BinWidth);
		}
		return Integral;
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
	char PMT_R11410::OnePhoton (Double_t time, PulseArray& electrons, bool fDebug) {
		Pulse OnePulse;     // Temporary variable for saving each SPE
		Double_t RND = 0;   // Random in range 0..1 defining interact. type
			// We can divide range 0..1 into 5 bands:
			// 0                    .. fProb_C1            - it was SPE from PC
			// fProb_C1             .. fProb_C             - it was DPE from PC
			// fProb_C              .. fProb_C + fProb_1d1 - it was SPE from 1dyn
			// fProb_C + fProb_1d1  .. fProb_C + fProb_1d  - it was DPE from 1dyn
			// fProb_C + fProb_1d   .. 1                   - it was no interaction
		fRND.SetSeed(0);
		RND = fRND.Rndm();
		char InteractType = 0; //0-no int, 1 - 1phe in PC, 2 - 2phe in PC, -1 - 1phe in 1d, -2 - 2phe in 1d

		if (fDebug) cout << "At " << time/ns << " ns ";

		// If photon didn't interact with PMT, nothing to do 
		if (RND > (fProb_C + fProb_1d)) {
			if (fDebug) cout << "photon didn't interact with PMT" << endl;
			InteractType = 0;
		}

		// If photon created 1 phe in PC
		else if (RND < fProb_C1) {
			if (fDebug) cout << "photon created 1 phe in PC" << endl;
			InteractType = 1;
			OnePulse.fAmpl = fRND.Gaus (fAmpl_mean, fAmpl_sigma);
			OnePulse.fTime = fRND.Gaus (fTOFe_mean, fTOFe_sigma) + time;
			if (fDebug) cout << "with amplitude = " << OnePulse.fAmpl << " and time from '0' of OutWave =" << OnePulse.fTime/ns << " ns" << endl;
			electrons.push_back (OnePulse);
		}

		// If photon created 2 phe in PC
		else if (RND < fProb_C) {
			if (fDebug) cout << "photon created 2 phe in PC" << endl;
			InteractType = 2;
			for (int i = 0; i < 2; i++) {
				OnePulse.fAmpl = fRND.Gaus (fAmpl_mean, fAmpl_sigma);
				OnePulse.fTime = fRND.Gaus (fTOFe_mean, fTOFe_sigma) + time;
				if (fDebug) cout << "with amplitude = " << OnePulse.fAmpl << " and time from '0' of OutWave =" << OnePulse.fTime/ns << " ns" << endl;
				electrons.push_back (OnePulse);
			}           
		}

		// If photon created 1 phe in 1dyn
		else if (RND < (fProb_C + fProb_1d1)) {
			if (fDebug) cout << "photon created 1 phe in 1dyn" << endl;
			InteractType = -1;
			OnePulse.fAmpl = fRND.Gaus (fAmpl_mean/fGain_PC_1d, fAmpl_sigma);
			OnePulse.fTime = fRND.Gaus (fTOFe_1d_mean, fTOFe_1d_sigma) + time;
			if (fDebug) cout << "with amplitude = " << OnePulse.fAmpl << " and time from '0' of OutWave =" << OnePulse.fTime/ns << " ns" << endl;
			electrons.push_back (OnePulse);
		}

		// Else photon created 2 phe in 1dyn
		else {
			if (fDebug) cout << "photon created 2 phe in 1dyn" << endl;
			InteractType = -2;
			for (int i = 0; i < 2; i++) {
				OnePulse.fAmpl = fRND.Gaus (fAmpl_mean/fGain_PC_1d, fAmpl_sigma);
				OnePulse.fTime = fRND.Gaus (fTOFe_1d_mean, fTOFe_1d_sigma) + time;
				if (fDebug) cout << "with amplitude = " << OnePulse.fAmpl << " and time from '0' of OutWave =" << OnePulse.fTime/ns << " ns" << endl;
				electrons.push_back (OnePulse);
			}
		}

		return (InteractType);
	}
	
	void PMT_R11410::Clear(Option_t *option) { 
		cout << "Function 'Clear' is not yet implemented" << endl;
	}

	void PMT_R11410::Print(Option_t *option) const {
		cout << endl << "Printing PMT parameters..." << endl;
		cout << "Defined parameters:" << endl;
		cout << "  SPE shape type   =  ";
		switch (fMode) {
			case 0 :
				cout << "None" << endl;
				break;
			case 1 :
				cout << "TF1" << endl;
				break;
			case 2 :
				cout << "TSpline" << endl;
				break;
		}
		cout <<    "  SPE Xmin         =  " << GetXmin()/ns << " ns";
		cout << " \t//SPE shape Xmin"       << endl;
		cout <<    "  SPE Xmax         =  " << GetXmax()/ns << " ns";
		cout << " \t//SPE shape Xmax"       << endl;
		cout <<    "  SPE Area mean    =  " << fArea_mean/(mV*ns) << " mV*ns";
		cout << " \t//SPE Area, mean"       << endl;
		cout <<    "  SPE Area sigma   =  " << fArea_sigma/(mV*ns) << " mV*ns";
		cout << " \t//SPE Area, sigma"      << endl;
		cout <<    "  QE of PC         =  " << fQE*100 << " %";
		cout << " \t//Quantum Efficiency of Photocathode" << endl;
		cout <<    "  QE of 1d         =  " << fQE_1d*100 << " %";
		cout << " \t//Quantum Efficiency of 1st dynode" << endl;
		cout <<    "  DPE_PC           =  " << fDPE_PC*100 << " %";
		cout << " \t//Double Photoelectron Emission probability for PC" << endl;
		cout <<    "  DPE_1d           =  " << fDPE_1d*100 << " %";
		cout << " \t//Double Photoelectron Emission probability for 1st dynode" << endl;
		cout <<    "  GF of 1d         =  " << fGF_1d*100 << " %";
		cout << " \t//Geometric factor (or probability for a random photon from PC to hit 1st dynode)" << endl;
		cout <<    "  PC->1d gain      =  " << fGain_PC_1d;
		cout << " \t//Amplification in the Photocathode - 1st dynode gap" << endl;
		cout <<    "  DCR              =  " << fDCR << " Hz";
		cout << " \t//Dark count rate"      << endl;
		cout <<    "  TOF PC->An mean  =  " << fTOFe_mean/ns << " ns";
		cout << " \t//Time of Flight from Photocathode to anode, mean" << endl;
		cout <<    "  TOF PC->An sigma =  " << fTOFe_sigma/ns << " ns";
		cout << " \t//Time of Flight from Photocathode to anode, sigma" << endl;
		cout <<    "  TOF PC->1d mean  =  " << fTOFe_PC_1d/ns << " ns";
		cout << " \t//Time of Flight from Photocathode to 1st dynode, mean" << endl;
		cout <<    "  AP cont          =  " << fAP_cont;
		cout << " \t//Afterpulses ?"        << endl;
		cout <<    "  AP peak          =  " << fAP_peak;
		cout << " \t//Afterpulses ?"        << endl;
		
		cout << "Calculated parameters:"        << endl;
		cout << std::setiosflags(std::ios::fixed) << std::setprecision(2);
		cout <<    "  Total PC probability =  " << fProb_C*100 << " %";
		cout << " \t//Total interaction probability on Photocathode" << endl;
		cout <<    "  1phe PC probability  =  " << fProb_C1*100 << " %";
		cout << " \t//1phe interaction probability on Photocathode" << endl;
		cout <<    "  2phe PC probability  =  " << fProb_C2*100 << " %";
		cout << " \t//2phe interaction probability on Photocathode" << endl;
		cout <<    "  Total 1d probability =  " << fProb_1d*100 << " %";
		cout << " \t//Total interaction probability on 1st dynode" << endl;
		cout <<    "  1phe 1d probability  =  " << fProb_1d1*100 << " %";
		cout << " \t//Total interaction probability on 1st dynode" << endl;
		cout <<    "  2phe 1d probability  =  " << fProb_1d2*100 << " %";
		cout << " \t//Total interaction probability on 1st dynode" << endl;
		cout << std::resetiosflags(std::ios::fixed) << std::setprecision(6);
		cout <<    "  SPE 1d Area mean     =  " << fArea_1d_mean/(mV*ns) << " mV*ns";
		cout << " \t//1st dynode SPE Area, mean" << endl;
		cout <<    "  SPE 1d Area sigma    =  " << fArea_1d_sigma/(mV*ns) <<" mV*ns";
		cout << " \t//1st dynode SPE Area, sigma" << endl;
		cout <<    "  TOF 1d->An mean      =  " << fTOFe_1d_mean/ns << " ns";
		cout << " \t//Time of Flight from 1st dynode to anode, mean" << endl;
		cout <<    "  TOF 1d->An sigma     =  " << fTOFe_1d_sigma/ns << " ns";
		cout << " \t//Time of Flight from 1st dynode to anode, sigma" << endl;
		if (fMode == 1) {
			cout <<    "  SPE Shape Area       =  " << fShape.func->Integral(GetXmin (), GetXmax ())/(ns) << " arb.un.*ns";
			cout << " \t//SPE shape Area, integrated from Xmin to Xmax by ROOT standard method" << endl;
			cout <<    "  SPE Shape MinValue   =  " << GetYmin() << " arb.un.";
			cout << " \t//Minimum value of SPE shape in range Xmin .. Xmax" << endl;
			cout <<    "  SPE Shape MaxValue   =  " << GetYmax() << " arb.un.";
			cout << " \t//Maximum value of SPE shape in range Xmin .. Xmax" << endl;
		}
	}
}