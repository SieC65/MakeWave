#ifndef PMT_R11410_HH
#define PMT_R11410_HH
#include <TF1.h>
#include <TSpline.h>
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

namespace CLHEP 
{
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
			
			// Structure for pulse parameters and vector of pulses
			struct Pulse {
				Double_t fAmpl;
				Double_t fTime;
			};
			typedef std::vector<struct Pulse> PulseArray;

			// The enumerator lists possible types of SPE shape
			enum {
				kModeNone,
				kModeF1,
				kModeSpline
			} fMode;

			// The union contain pointer to SPE shape in form of TF1 or TSpline
			union SPEShape {
				TF1     *func;
				TSpline *spline;
			} fShape;

			// SETTERS
			virtual void SetDefaults() = 0; // Default values for all parameters
			virtual void SetParams (double QE      , double Area_mean, 
			                        double DCR = 0., double AP_cont = 0.) = 0; // Main control parameters
			virtual void SetShape (TF1     *Shape) = 0; // SPE shape in function form
			virtual void SetShape (TSpline *Shape) = 0; // SPE shape in spline form

			// GETTERS
			virtual Double_t GetXmax()         const = 0; // Right point of the domain of definition
			virtual Double_t GetXmin()         const = 0; // Left point of the domain of definition
			virtual Double_t GetYmax()         const = 0; // Maximum value of SPE Shape
			virtual Double_t GetYmin()         const = 0; // Minimum value of SPE Shape
			virtual Double_t GetArea()         const = 0; // Pulse area of SPE
			virtual Double_t Eval (Double_t t) const = 0; // Value of SPE Shape at time t
			virtual Double_t GetShapeArea()    const = 0; // Pulse area of SPE Shape
			virtual Double_t GetAmpl()         const = 0;
			virtual Double_t GetAmpl_Sigma()   const = 0;

			// OUTPUT
			virtual void Print             (Option_t *option="") const { ; } // Print all PMT parameters
			virtual void PrintUsrDefParams (Option_t *option="") const { ; } // Print user-defined PMT parameters
			virtual void PrintCalcParams   (Option_t *option="") const { ; } // Print calculated PMT parameters
			virtual void PrintProbs        (Option_t *option="") const { ; } // Print interactions probabilities only
			virtual void DrawShape (char* title) const { ; } // Draw SPE shape
			
			// ACTIONS
			virtual int  Begin     (PulseArray &electrons) { return(0); }
			virtual Char_t OnePhoton (Double_t *time, PulseArray &Pulse, bool fDebug = false) = 0; // Convert photons to pulses
			virtual int  End       (PulseArray &electrons) { return(0); }
			virtual void Clear     (Option_t *option="") { ; }
			virtual void GenDCR    (Double_t begintime, Double_t endtime, PulseArray& DarkPulse) { ; } // Generate pulses for dark counts

		private:
		
			// Prevent simple copy
			PMT(const PMT &r);
			PMT & operator=(const PMT &r);
	};

	// Hamamatsu R11410-20 PMT
	class PMT_R11410 : public PMT
	{
			friend void test (PMT_R11410 &);
			
		public:
		
			 PMT_R11410();
			~PMT_R11410() {;}
			virtual TObject* Clone (const char *newname="") const { return 0;}
			
			// SETTERS
			void SetDefaults();
			void SetParams (double QE = 0.3, double Area_mean = 10*mV*ns, 
			                double DCR = 0., double AP_cont = 0.); // Main control parameters
			void SetShape (TF1     *Shape);
			void SetShape (TSpline *Shape);
			// Set advanced parameters
			void CalculateParams(); // Calculate auxiliary parameters
			void SetDPE_PC     (Double_t DPE_PC     = 0.225  );
			void SetDPE_1d     (Double_t DPE_1d     = 0      );
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
			Double_t GetXmax()        const;
			Double_t GetXmin()        const;
			Double_t GetYmax()        const;
			Double_t GetYmin()        const;
			Double_t GetArea()        const {return fArea_mean;}
			Double_t GetArea_Sigma()  const {return fArea_sigma;}
			Double_t GetShapeArea()   const {return fShapeArea;}
			Double_t GetAmpl()        const {return fAmpl_mean;}
			Double_t GetAmpl_Sigma()  const {return fAmpl_sigma;}
			Double_t Eval(Double_t t) const;
			// Get independent PMT parameters
			Double_t GetQE()          const {return fQE;}
			Double_t GetDPE_PC()      const {return fDPE_PC;}
			Double_t GetDPE_1d()      const {return fDPE_1d;}
			Double_t GetTOFe()        const {return fTOFe_mean;}
			Double_t GetTOFe_Sigma()  const {return fTOFe_sigma;}
			Double_t GetQE_1d()       const {return fQE_1d;}
			Double_t GetGain_PC_1d()  const {return fGain_PC_1d;}
			Double_t GetGF_1d()       const {return fGF_1d;}
			Double_t GetTOFe_PC_1d()  const {return fTOFe_PC_1d;}
			Double_t GetDCR()         const {return fDCR;}
			Double_t GetAP_cont()     const {return fAP_cont;}
			Double_t GetAP_peak()     const {return fAP_peak;}

			// ACTIONS
			int  Begin     (PulseArray &electrons);
			int  End       (PulseArray &electrons);
			Char_t OnePhoton (Double_t *time, PulseArray &electrons, bool fDebug=true);
			void GenDCR    (Double_t begintime, Double_t endtime, PulseArray& electrons);
			void Clear     (Option_t *option="");

			// OUTPUT INFO
			void Print             (Option_t *option="") const;
			void PrintUsrDefParams (Option_t *option="") const;
			void PrintCalcParams   (Option_t *option="") const;
			void PrintProbs        (Option_t *option="") const;
			void DrawShape (const char* title="SPE Shape") const;
			
		private:
		
			// User-defined PMT parameters
			Double_t fQE;         // Quantum Efficiency(standard, i.e. with many phe)
			Double_t fDPE_PC;     // Double Photoelectron Emission probability for PC
			Double_t fDPE_1d;     // Double Photoelectron Emission probability for 1dyn
			Double_t fArea_mean;  // SPE pulse area
			Double_t fArea_sigma; // Sigma for spreading SPE area by gauss
			Double_t fTOFe_mean;  // Time of Flight from photocathode to anode
			Double_t fTOFe_sigma; // Sigma for spreading ToF PC-anode by gauss
			Double_t fQE_1d;      // Quantum Efficiency for 1dyn
			Double_t fGain_PC_1d; // Amplification on first gap(cathode-1dynode), ratio
			Double_t fGF_1d;      // Average geom. prob. for a random photon from cathode to hit 1dyn
			Double_t fTOFe_PC_1d; // Time from cathode to 1dyn
			Double_t fDCR;        // Dark count rate, counts per nanosecond
			Double_t fAP_cont;    // Afterpulsing probability (continuum)
			Double_t fAP_peak;    // Afterpulsing probability (peak)

			// Auxiliary, used in calculations
			Double_t fShapeArea;       // SPE shape area
			Double_t fAmpl_mean;       // SPE Amplitude (nondimensional, relative to SPE Shape)
			Double_t fAmpl_sigma;      // Sigma of SPE amplitude for spreading by gauss
			Double_t fProb_C;          // Prob. of int. in PC
			Double_t fProb_C1;         // Prob. of 1phe in PC
			Double_t fProb_C2;         // Prob. of 2phe in PC
			Double_t fProb_1d;         // Full prob. of int. in 1dyn
			Double_t fProb_1d1;        // Full prob. of 1phe in 1dyn
			Double_t fProb_1d2;        // Full prob. of 2phe in 1dyn
			Double_t fArea_1d_mean;    // SPE area from 1d
			Double_t fArea_1d_sigma;   // = Area_sigma / Gain_PC_1d
			Double_t fTOFe_1d_mean;    // Time of Flight e- from 1dyn to anode
			Double_t fTOFe_1d_sigma;   // 0.5*TOF_sigma    
			TRandom3 fRND;             // Object for random calculations

			bool fDebug; //some extended info (just for debug)

			// ACTIONS
			// Calculate integral of TSpline object in range (xmin,xmax) by rectangles method
			Double_t SplineIntegral (TSpline* spline, Double_t xmin, Double_t xmax, Int_t nbins);
	};
  
} // namespace

#endif // PMT_R11410_HH
