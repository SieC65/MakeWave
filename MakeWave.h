#ifndef MakeWave_H
#define MakeWave_H
#include <vector>
#include <Rtypes.h>
#include <TF1.h>
#include "SystemOfUnits.h"

using namespace std;
using namespace CLHEP;

class MakeWave
{
	public:
		MakeWave ();
		void SetSPE (TF1 *SPE);		//Set SPE function
		void SetParams (Int_t Num, Double_t Delay, Double_t Period, Double_t Gain);	//Set parameters
		void SetTimeSeq (vector <double> *Tseq);	//Set sequence of SPE's arrival times
		void CreateOutWave ();	//Create OutWave with unit of ADC as parameter
		void PrintOutWave ();	//print all values of output signal
		void Draw ();			//Draw OutWave values from abs time
		vector <double> OutWave;	//Output signal
	private:
		vector <double> *ftimeseq;
		Double_t fPeriod;
		Double_t fDelay;
		Double_t fGain;
		Int_t fNum;
		TF1 *fSPE;
};

#endif // MakeWave_H
