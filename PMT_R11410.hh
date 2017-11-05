#ifndef PMT_R11410_HH
#define PMT_R11410_HH

#include <TF1.h>
#include <TSpline.h>
#include <TString.h>
#include <TRandom3.h>
#include "SystemOfUnits.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RED::PMT                                                             //
//                                                                      //
// A virtual class for PMT response to the arriving photons.            //
// It contains a pulse shape function for a Single-Photo-Electron       //
// response and number of pulses function that estimates amount and     //
// parameters of generated SPEs.                                        //
//                                                                      //
// Method OnePhoton() converts a single photon hit into a number of     //
// photoelectrons together with randomized amplitude and delay time.    //
// It can also take into account after-pulses.                          //
// Two additional methods Begin() and End() help to account for         //
// additional effects like dark count and saturation.                   //
//                                                                      //
// Method Eval() provides pulse shape function evaluation.              //
// The function must be finite itself and must have a finite interval   //
// as a domain of definition, yet it is recommended to be smooth.       //
// The amplitude and delay provided by OnePhoton() method can be used   //
// with this function in the form fAmplitude*Eval(time-fDelay),         //
// resulting in waveform signal. Time and amplitude here use units      //
// from CLHEP/SystemOfUnits.h or compatible.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

namespace CLHEP {
	static const double mV = 1.e-3*CLHEP::volt;
}
using CLHEP::ns;
using CLHEP::mV;

namespace RED
{
	class PMT : public TObject
	{
	public:
		PMT() { ; }
		virtual ~PMT() { ; }
		//virtual TObject* Clone(const char *newname="") const = 0;

		// OUTPUT INFO
		virtual void Print (Option_t *option="") const { ; } // Print all PMT parameters
		virtual void DrawShape (char* title)     const { ; } // Draw SPE shape

		// SETTERS
		virtual void SetDefaults() = 0; // Default values for all parameters
		virtual void SetParams (double QE      , double Area_mean, 
		                        double DCR = 0., double AP_cont = 0.) = 0; // Main control parameters
		virtual void SetShape (TF1     *Shape) = 0; // SPE shape in function form
		virtual void SetShape (TSpline *Shape) = 0; // SPE shape in spline form

		// GETTERS
		// Get SPE parameters
		virtual Double_t GetXmax()         const = 0; // Right point of the domain of definition
		virtual Double_t GetXmin()         const = 0; // Left point of the domain of definition
		virtual Double_t GetYmax()         const = 0; // Maximum value of SPE shape
		virtual Double_t GetYmin()         const = 0; // Minimum value of SPE shape
		virtual Double_t GetArea()         const = 0; // Pulse area of SPE shape
		virtual Double_t Eval (Double_t t) const = 0; // Value of SPE shape

		// Structure for pulse parameters and vector of pulses
		struct Pulse {
			Double_t fAmpl;
			Double_t fTime;
		};
		typedef std::vector<struct Pulse> PulseArray;

		// MAIN FUNCTIONS
		virtual int Begin (PulseArray &electrons) { return(0); }
		virtual int OnePhoton (Double_t time, PulseArray &electrons,
		                       bool fDebug = false) = 0; // Convert photon to pulses
		virtual int End (PulseArray &electrons) { return(0); }
		virtual void Clear (Option_t *option="") { ; }

		// Define SPE shape type
		enum {
			kModeNone,
			kModeF1,
			kModeSpline
		} fMode;

		// Allocate SPE shape type
		union SPEShape {
			TF1     *func;
			TSpline *spline;
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
		virtual TObject* Clone (const char *newname="") const { return 0;}

		// OUTPUT INFO
		void Print (Option_t *option="") const;
		void DrawShape (const char* title="SPE Shape");

		// SETTERS
		void SetDefaults();
		void SetParams (double QE = 0.3, double Area_mean = 10*mV*ns, 
		                double DCR = 0., double AP_cont = 0.); // Main control parameters
		void SetShape (TF1     *Shape);
		void SetShape (TSpline *Shape);
		// Set advanced parameters
		void CalculateParams(); // Calculate auxiliary parameters
		void SetDPE        (Double_t DPE        = 0.225  );
		void SetQE_1d      (Double_t QE_1d      = 0.105  );
		void SetGain_PC_1d (Double_t Gain_PC_1d = 13     );
		void SetGF_1d      (Double_t GF_1d      = 0.1    );
		void SetArea_sigma (Double_t Area_sigma = 1*mV*ns);
		void SetTOFe_PC_1d (Double_t TOFe_PC_1d = 6*ns   );
		void SetTOFe_mean  (Double_t TOFe_mean  = 30*ns  );
		void SetTOFe_sigma (Double_t TOFe_sigma = 3*ns   );
		void SetAP_peak    (Double_t AP_peak    = 0      );


		// GETTERS
		// Get SPE parameters
		Double_t GetXmax()            const;
		Double_t GetXmin()            const;
		Double_t GetYmax()            const;
		Double_t GetYmin()            const;
		Double_t GetArea()            const {return fArea_mean;}
		Double_t GetArea_Sigma()      const {return fArea_sigma;}
		Double_t GetAmpl()            const {return fAmpl_mean;}
		Double_t GetAmpl_Sigma()      const {return fAmpl_sigma;}
		Double_t Eval(Double_t t)     const;
		// Get independent PMT parameters
		Double_t GetQE()              const {return fQE;}
		Double_t GetDPE()             const {return fDPE;}
		Double_t GetTOFe()            const {return fTOFe_mean;}
		Double_t GetTOFe_Sigma()      const {return fTOFe_sigma;}
		Double_t GetQE_1d()           const {return fQE_1d;}
		Double_t GetGain_PC_1d()      const {return fGain_PC_1d;}
		Double_t GetGF_1d()           const {return fGF_1d;}
		Double_t GetTOFe_PC_1d()      const {return fTOFe_PC_1d;}
		Double_t GetDCR()             const {return fDCR;}
		Double_t GetAP_cont()         const {return fAP_cont;}
		Double_t GetAP_peak()         const {return fAP_peak;}

		// MAIN FUNCTIONS
		int Begin (PulseArray &electrons);
		int OnePhoton (Double_t time, PulseArray &electrons, bool fDebug=true);
		int End   (PulseArray &electrons);
		void Clear (Option_t *option="");

	private:
		// Native parameters
		// Area and TOF comes from making several photoelectrons

		// Independent PMT parameters
		Double_t fQE;         // Quantum Efficiency(standard, i.e. with many phe), probability
		Double_t fDPE;        // Double Photoelectron Emission P(2phe)/(P(2phe)+P(1phe))
		Double_t fArea_mean;  // SPE pulse area
		Double_t fArea_sigma; //
		Double_t fTOFe_mean;  // Time of Flight from cathode to anode
		Double_t fTOFe_sigma; //
		Double_t fQE_1d;      // Quantum Efficiency for 1dyn(only 1phe), probability
		Double_t fGain_PC_1d; // Amplification on first gap(cathode-1dynode), ratio
		Double_t fGF_1d;      // Average geom. prob. for a random photon from cathode to hit 1st dynode
		Double_t fTOFe_PC_1d; // Time from cathode to 1st dynode
		Double_t fDCR;        // Dark count rate, Hz
		Double_t fAP_cont;    //
		Double_t fAP_peak;    // Afterpulsing probability, for continuum and peak

		// Auxiliary, actually used in calculation
		Double_t fAmpl_mean;       // SPE Amplitude
		Double_t fAmpl_sigma;      //
		Double_t fProb_C;          // Prob. of int. in PC
		Double_t fProb_C1;         // Prob. of 1phe in PC
		Double_t fProb_C2;         // Prob. of 2phe in PC
		Double_t fProb_1d;         // Full prob. of int. in 1d
		Double_t fArea_1d_mean;    // SPE area from 1d
		Double_t fArea_1d_sigma;   // = Area_sigma / Gain_PC_1d
		Double_t fTOFe_1d_mean;    // Time of Flight e- from 1dyn to anode
		Double_t fTOFe_1d_sigma;   // 0.5*TOF_sigma    
		TRandom3 fRND;             //

		bool fDebug; //some extended info (just for debug)
		Double_t SplineIntegral (TSpline* spline, Double_t xmin, Double_t xmax, Int_t nbins);
	};
  
} // namespace

#endif // PMT_R11410_HH