#include <iostream>
#include <vector>
#include <cmath>

#include "Pythia8/Pythia.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"


const double gSMEARING_FACTOR = 0.001;

struct particleData {
	
	double E;
	double px;
	double py;
	double pz;
	int    nEvent;
};

int main(int argc, char* argv[]) {

	const int nEvents = atoi(argv[1]);
	
	TFile *outfile = new TFile("trees.root", "RECREATE");

	TRandom3 *rand = new TRandom3();

	particleData particle;

	TTree *gammasTree    = new TTree("GammasTree", "GammasTree");
	TTree *piPlusesTree  = new TTree("PiPlusesTree", "PiPlusesTree");
	TTree *piMinusesTree = new TTree("PiMinusesTree", "PiMinusesTree");

	gammasTree   ->Branch("ParticleData", &particle, "E/D:px/D:py/D:pz/D:nEvent/I");
	piPlusesTree ->Branch("ParticleData", &particle, "E/D:px/D:py/D:pz/D:nEvent/I");
	piMinusesTree->Branch("ParticleData", &particle, "E/D:px/D:py/D:pz/D:nEvent/I");

	Pythia8::Pythia pythia;
	//pythia.readString("HardQCD:all = on");
	//pythia.readString("PhaseSpace:pTHatMin = 20.");
	pythia.readString("SoftQCD:all = on");
	pythia.readString("Beams:eCM = 200.");

	pythia.init();


	for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

  		if (!pythia.next()) continue;

    		for (int i = 0; i < pythia.event.size(); ++i) {
			
			if (pythia.event[i].isFinal() && pythia.event[i].id() == 22) {

				particle.E  = rand->Gaus(pythia.event[i].e(),  gSMEARING_FACTOR);		
				particle.px = rand->Gaus(pythia.event[i].px(), gSMEARING_FACTOR);		
				particle.py = rand->Gaus(pythia.event[i].py(), gSMEARING_FACTOR);		
				particle.pz = rand->Gaus(pythia.event[i].pz(), gSMEARING_FACTOR);		
				particle.nEvent = iEvent;

				gammasTree   ->Fill();
			}

			if (pythia.event[i].isFinal() && pythia.event[i].id() == 211) {

				particle.E  = rand->Gaus(pythia.event[i].e(),  gSMEARING_FACTOR);		
				particle.px = rand->Gaus(pythia.event[i].px(), gSMEARING_FACTOR);		
				particle.py = rand->Gaus(pythia.event[i].py(), gSMEARING_FACTOR);		
				particle.pz = rand->Gaus(pythia.event[i].pz(), gSMEARING_FACTOR);		
				particle.nEvent = iEvent;

				piPlusesTree->Fill();
			}

			if (pythia.event[i].isFinal() && pythia.event[i].id() == -211) {

				particle.E  = rand->Gaus(pythia.event[i].e(),  gSMEARING_FACTOR);		
				particle.px = rand->Gaus(pythia.event[i].px(), gSMEARING_FACTOR);		
				particle.py = rand->Gaus(pythia.event[i].py(), gSMEARING_FACTOR);		
				particle.pz = rand->Gaus(pythia.event[i].pz(), gSMEARING_FACTOR);		
				particle.nEvent = iEvent;

				piMinusesTree->Fill();
			}
			
		}

	}
	
	gammasTree   ->Write();
	piMinusesTree->Write();
	piPlusesTree ->Write();

	outfile->Close();

 	pythia.stat();

	return 0;
}


