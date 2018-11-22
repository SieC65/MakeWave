#include <algorithm>
#include <vector>

#include <Rtypes.h>

#include "MakeWave.h"

class MakeTest
{
	public:

		MakeTest();

		//SETTERS
		void SetPhotonsTimes (Double_t *TimesArr); // Set vector of photon arrival times

		//GETTERS
		vector <Double_t> GetMeanOW  () {return *fMeanOutWave;}
		vector <Double_t> GetTimesVec() {return *fTimesVec;}

		//ACTIONS
		void AverageOW (Int_t DebugN, MakeWave* MakeWaveObj); // Create average OutWave for DebugN OutWaves with the same SPE arrival times
		void DrawOW (MakeWave* MakeWaveObj);
		void RandPhotonTimes (Int_t number, Double_t leftEdge, Double_t rightEdge); // Set random photon arrival times
		void PrintTimesVec ();
		void DrawShape (RED::PMT *pmt, Char_t const *title);

	private:

		vector <Double_t> *fTimesVec; // vector of arrival times
		vector <Double_t> *fMeanOutWave; // outwave vector
};