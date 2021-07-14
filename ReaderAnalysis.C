#include <cmath>

struct particleData {
	double E;
	double px;
	double py;
	double pz;
	int    nEvent;
};

int ReaderAnalysis() {

	double M;

	int nEvents;

	bool check = true;

	particleData gamma, gammaAux, piPlus, piMinus;		

	std::vector<particleData> event;

	TFile *treesFile = new TFile("trees.root",   "READ");
	TFile *outfile   = new TFile("outfile.root", "RECREATE");

	TTreeReader *gammaReader   = new TTreeReader("GammasTree",      treesFile);
	TTreeReader *piPlusReader  = new TTreeReader("PiPlusesTree",  treesFile);
	TTreeReader *piMinusReader = new TTreeReader("PiMinusesTree", treesFile);

	TTreeReaderValue<double> gammaE      (*gammaReader, "ParticleData.E");
	TTreeReaderValue<double> gammaPx     (*gammaReader, "ParticleData.px");
	TTreeReaderValue<double> gammaPy     (*gammaReader, "ParticleData.py");
	TTreeReaderValue<double> gammaPz     (*gammaReader, "ParticleData.pz");
	TTreeReaderValue<int>    gammaNevent (*gammaReader, "ParticleData.nEvent");

	nEvents = treesFile->Get<TTree>("GammasTree")->GetMaximum("ParticleData.nEvent") + 1;

	TH1D *piEtaHist = new TH1D("PiEtaHist", "PiEtaHist", 500,0,800);

	gammaReader->Next();

	for (int iEvent = 0; iEvent < nEvents; iEvent++) {
		
		while(1) {

			gamma.E      = *gammaE;
			gamma.px     = *gammaPx;
			gamma.py     = *gammaPy;
			gamma.pz     = *gammaPz;
			gamma.nEvent = *gammaNevent;
		
			event.push_back(gamma);

			check = gammaReader->Next();

			std::cout << iEvent << "   " << gamma.E << "   " << check << std::endl;

			if (*gammaNevent != iEvent) break;
			if ( !check )               break;
		}

		for (int iGamma = 0; iGamma < event.size(); iGamma++) {
			for (int jGamma = iGamma + 1; jGamma < event.size(); jGamma++) {

				M = 1000 * sqrt( (event[iGamma].E  + event[jGamma].E)  * (event[jGamma].E  + event[jGamma].E) 
						-(event[iGamma].px + event[jGamma].px) * (event[jGamma].px + event[jGamma].px) 
						-(event[iGamma].py + event[jGamma].py) * (event[jGamma].py + event[jGamma].py)
						-(event[iGamma].pz + event[jGamma].pz) * (event[jGamma].pz + event[jGamma].pz) );

				piEtaHist->Fill(M);
			}
		}

		event.clear();
	}


	piEtaHist->Write();

	return 0;
}
