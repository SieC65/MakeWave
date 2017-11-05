#ifndef PMT_R11410_HH
#define PMT_R11410_HH

#include <TF1.h>
#include <TSpline.h>
#include <TString.h>
#include <TRandom3.h>

//////////////////////////////////////////////////////////////////////////
//																		//
// RED::PMT																//
//																		//
// A virtual class for PMT response to the arriving photons.			//
// It contains a pulse shape function for a Single-Photo-Electron	 	//
// response and number of pulses function that estimates amount and		//
// parameters of generated SPEs.										//
//																		//
// Method OnePhoton() converts a single photon hit into a number of		//
// photoelectrons together with randomized amplitude and delay time.	//
// It can also take into account after-pulses.							//
// Two additional methods Begin() and End() help to account for			//
// additional effects like dark count and saturation.					//
//																		//
// Method Eval() provides pulse shape function evaluation.				//
// The function must be finite itself and must have a finite interval	//
// as a domain of definition, yet it is recommended to be smooth.		//
// The amplitude and delay provided by OnePhoton() method can be used	//
// with this function in the form fAmplitude*Eval(time-fDelay),			//
// resulting in waveform signal. Time and amplitude here use units		//
// from CLHEP/SystemOfUnits.h or compatible.							//
//																		//
//////////////////////////////////////////////////////////////////////////

namespace RED
{

  class PMT : public TObject
  {
  public:
	PMT() { ; }
	virtual ~PMT() { ; }
	//virtual TObject* Clone(const char *newname="") const = 0;

	// OUTPUT INFO
	virtual void Print(Option_t *option="") const { ; } // Print all PMT parameters
	virtual void DrawShape(char* title) const { ; } // Draw SPE shape

	// SETTERS
	virtual void SetDefaults() = 0; // Some default values for all parameters
	virtual void SetParameters(	double_t QE, double_t DPE, double_t QE_1d, double_t Gain_PC_1d, 
								double_t GF_1d, double_t Area_mean, double_t Area_sigma, 
								double_t TOFe_PC_1d, double_t TOFe_mean, double_t TOFe_sigma, 
								double_t DCR, double_t AP_cont, double_t AP_peak) = 0; // Main control parameters
	virtual void SetShape(TF1 *Shape) = 0; 		// Shape in function form
	virtual void SetShape(TSpline *Shape) = 0; 	// Shape in spline form

	// GETTERS
	// Get SPE parameters
	virtual Double_t GetXmax() 			const = 0; // Right point of the domain of definition
	virtual Double_t GetXmin() 			const = 0; // Left point of the domain of definition
	virtual Double_t GetMaximum() 		const = 0; // The maximum value of the function
	virtual Double_t GetMinimum() 		const = 0; // The mimimum value of the function
	virtual Double_t Eval(Double_t t) 	const = 0; // SPE pulse shape function
	// Get independent PMT parameters
	virtual Double_t GetQE() 			const = 0; // Quantum Efficiency(standard, i.e. with many phe), probability
	virtual Double_t GetArea() 			const = 0; // The pulse area in the function measured as for SPE
	virtual Double_t GetArea_Sigma() 	const = 0;
	virtual Double_t GetAmplitude() 	const = 0; // SPE Amplitude
	virtual Double_t GetAmplitude_Sigma() const = 0;
	virtual Double_t GetTOFe() 			const = 0; // Time of Flight from cathode to anode, s
	virtual Double_t GetTOFe_Sigma() 	const = 0;
	virtual Double_t GetDPE() 			const = 0; // Double Photoelectron Emission P(2phe)/(P(2phe)+P(1phe))
	virtual Double_t GetQE_1d() 		const = 0; // Quantum Efficiency for 1dyn(only 1phe), probability
	virtual Double_t GetGain_PC_1d() 	const = 0; // Amplification on first gap(cathode-1dynode), ratio
	virtual Double_t GetGF_1d() 		const = 0; // Average geometrical probability for a random photon from cathode to hit 1st dynode, probability
	virtual Double_t GetTOFe_PC_1d() 	const = 0; // Time from cathode to 1st dynode, s
	virtual Double_t GetDCR() 			const = 0; // Dark count rate, Hz
	virtual Double_t GetAP_cont() 		const = 0;
	virtual Double_t GetAP_peak() 		const = 0; // Afterpulsing probability, for continuum and peak
	// Get calculated PMT parameters
	virtual Double_t GetProb_C() 		const = 0; // Prob. of int. in PC
	virtual Double_t GetProb_C1() 		const = 0; // Prob. of 1phe in PC
	virtual Double_t GetProb_C2() 		const = 0; // Prob. of 2phe in PC
	virtual Double_t GetProb_1d() 		const = 0; // Full prob. of int. in 1d
	virtual Double_t GetArea_1d() 		const = 0; // SPE area from 1d
	virtual Double_t GetArea_1d_Sigma() const = 0;
	virtual Double_t GetTOFe_1d() 		const = 0; // Time of Flight e- from 1dyn to anode
	virtual Double_t GetTOFe_1d_Sigma() const = 0;
	
	// Structure for pulse parameters and vector of pulses
	struct Pulse {
		double_t fAmplitude;
		double_t fTime;
	};
	typedef std::vector<struct Pulse> PulseArray;
	
	// MAIN FUNCTIONS
	virtual int Begin(PulseArray &electrons) { return(0); }
	virtual int OnePhoton(double_t time, PulseArray &electrons, bool fDebug=false) = 0; // Convert photon to pulses
	virtual int End(PulseArray &electrons) { return(0); }
	virtual void Clear(Option_t *option="") { ; }
		
	// Define SPE shape type
	enum {
		kModeNone,
		kModeF1,
		kModeSpline
	} fMode;
	
	// Allocate SPE shape type
	union SPEShape {
		TF1 *f1;
		TSpline *f2;
	} fShape;

  private:
	// Prevent simple copy
	PMT(const PMT &r);
	PMT & operator=(const PMT &r);
  };

  // Hamamatsu R11410-20 PMT
  class PMT_R11410 : public PMT
  {
  public:
	PMT_R11410();
	~PMT_R11410() {;}
	virtual TObject* Clone(const char *newname="") const { return 0;}
	
	void Clear(Option_t *option="");
	void Print(Option_t *option="") const;
	void DrawShape(const char* title="SPE Shape");

	// SETTERS
	void SetDefaults();	// Set default parameters, such as QE, DPE, Areas, ToF, DCR
	void SetParameters(	double_t QE, double_t DPE, double_t QE_1d, double_t Gain_PC_1d, 
						double_t GF_1d, double_t Area_mean, double_t Area_sigma, 
						double_t TOFe_PC_1d, double_t TOFe_mean, double_t TOFe_sigma, 
						double_t DCR, double_t AP_cont, double_t AP_peak);
	void SetShape(TF1 *Shape);
	void SetShape(TSpline *Shape);
	void CalculateParams(); // Calculate auxiliary parameters

	// GETTERS
	// Get SPE parameters
	Double_t GetXmax() 			const;
	Double_t GetXmin() 			const;
	Double_t GetMaximum() 		const;
	Double_t GetMinimum() 		const;
	Double_t Eval(Double_t t) 	const;
	// Get independent PMT parameters
	Double_t GetQE() 			const {return fQE;}
	Double_t GetArea() 			const {return fArea_mean;}
	Double_t GetArea_Sigma() 	const {return fArea_sigma;}
	Double_t GetAmplitude() 	const {return fAmplitude_mean;}
	Double_t GetAmplitude_Sigma() const {return fAmplitude_sigma;}
	Double_t GetTOFe() 			const {return fTOFe_mean;}
	Double_t GetTOFe_Sigma() 	const {return fTOFe_sigma;}
	Double_t GetDPE() 			const {return fDPE;}
	Double_t GetQE_1d() 		const {return fQE_1d;}
	Double_t GetGain_PC_1d() 	const {return fGain_PC_1d;}
	Double_t GetGF_1d() 		const {return fGF_1d;}
	Double_t GetTOFe_PC_1d() 	const {return fTOFe_PC_1d;}
	Double_t GetDCR() 			const {return fDCR;}
	Double_t GetAP_cont() 		const {return fAP_cont;}
	Double_t GetAP_peak() 		const {return fAP_peak;}
	// Get calculated PMT parameters
	Double_t GetProb_C() 		const {return fProb_C;}
	Double_t GetProb_C1() 		const {return fProb_C1;}
	Double_t GetProb_C2() 		const {return fProb_C2;}
	Double_t GetProb_1d() 		const {return fProb_1d;}
	Double_t GetArea_1d() 		const {return fArea_1d_mean;}
	Double_t GetArea_1d_Sigma() const {return fArea_1d_sigma;}
	Double_t GetTOFe_1d() 		const {return fTOFe_1d_mean;}
	Double_t GetTOFe_1d_Sigma() const {return fTOFe_1d_sigma;}

	int Begin(PulseArray &electrons);
	int OnePhoton(double_t time, PulseArray &electrons, bool fDebug=true);
	int End(PulseArray &electrons);

  private:
	// Native parameters
	// Area and TOF comes from making several photoelectrons

	// Independent PMT parameters
	double_t fQE;				// Quantum Efficiency(standard, i.e. with many phe), probability
	double_t fDPE;				// Double Photoelectron Emission P(2phe)/(P(2phe)+P(1phe))
	double_t fArea_mean;		// SPE pulse area
	double_t fArea_sigma;		//
	double_t fTOFe_mean;		// Time of Flight from cathode to anode
	double_t fTOFe_sigma;		//
	double_t fQE_1d;			// Quantum Efficiency for 1dyn(only 1phe), probability
	double_t fGain_PC_1d;		// Amplification on first gap(cathode-1dynode), ratio
	double_t fGF_1d;			// Average geom. prob. for a random photon from cathode to hit 1st dynode
	double_t fTOFe_PC_1d;		// Time from cathode to 1st dynode
	double_t fDCR;				// Dark count rate, Hz
	double_t fAP_cont;			//
	double_t fAP_peak;			// Afterpulsing probability, for continuum and peak

	// Auxiliary, actually used in calculation
	double_t fAmplitude_mean;	// SPE Amplitude
	double_t fAmplitude_sigma;	//
	double_t fProb_C;			// Prob. of int. in PC
	double_t fProb_C1;			// Prob. of 1phe in PC
	double_t fProb_C2;			// Prob. of 2phe in PC
	double_t fProb_1d;			// Full prob. of int. in 1d
	double_t fArea_1d_mean;		// SPE area from 1d
	double_t fArea_1d_sigma;	// = Area_sigma / Gain_PC_1d
	double_t fTOFe_1d_mean;		// Time of Flight e- from 1dyn to anode
	double_t fTOFe_1d_sigma;	// 0.5*TOF_sigma	
	TRandom3 fRND;				//
	
	bool fDebug;	//some extended info (just for debug)
  };
  
} // namespace

#endif // PMT_R11410_HH