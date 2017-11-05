#ifndef PMT_H
#define PMT_H
#include <Rtypes.h>

class PMT
{
	public:
		PMT ();
		
		//Setters
		void SetPMTPar (Double_t QE, Double_t DPE, Double_t QE_1d, Double_t ArbAmpl_1d, Double_t GF_1d, Double_t AmplSigma, Double_t TOF_1d);	//Set PMT parameters
		
		//Getters
		Double_t GetPphe_PC ();
		Double_t GetPphe_1d ();
		Double_t GetArbAmpl_1d ();
		Double_t GetAmplSigma_1d ();
		Double_t GetDPE();
		Double_t GetTOF_1d();
	private:
		//PMT Parameters 
		//Set by user:
		Double_t fQE;		//Quantum Efficiency (full)
		Double_t fDPE;		//Double Photoelectron Emission P(2phe)/(P(2phe)+P(1phe))
		Double_t fQE_1d;	//Quantum Efficiency for 1dyn (only 1phe)
		Double_t fArbAmpl_1d;	//Ampl(1d)/Ampl(PC) for SPE
		Double_t fGF_1d;		//Arb. geometrical factor of 1dyn
		Double_t fTOF_1d;		//TOF for photon from PC to 1dyn
		//Calculate:
		Double_t fQE_1d_ratio;	//auxiliary quantity	
		Double_t fPphe_PC;		//Probability of photoeffect on PC
		Double_t fPphe_1d;		//Prob-ty of ph.effect on 1dyn for photons passed PC
		Double_t fAmplSigma_1d;	//Sigma for amplitude from 1dyn		
};

#endif // PMT_H
