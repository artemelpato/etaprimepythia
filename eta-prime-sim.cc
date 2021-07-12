#include <iostream>
#include <vector>
#include <cmath>

#include "Pythia8/Pythia.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"

const double gETA_MASS  = 547; //Mev
const double gETA_SIGMA = 15; //Mev
const double gSMEARING_FACTOR = 0.001;

struct fourVec {
	
	double E;
	double px;
	double py;
	double pz;
};

void FindEtaGammas( std::vector<fourVec>&,
        	    std::vector<fourVec>&,
            	    std::vector<fourVec>&,
         	    TH1D* );

void MixParticle( std::vector<fourVec>&,
                  std::vector<fourVec>&, 
        	  std::vector<fourVec>&, 
        	  std::vector<fourVec>&, 
        	  TH1D* );    

void MixParticleNew( std::vector<fourVec>&, 
		     std::vector<fourVec>&, 
		     std::vector<fourVec>&, 
		     TH1D* );    


int main(int argc, char* argv[]) {

	const int nEvents = atoi(argv[1]);
	
	TFile *outfile = new TFile("eta_sim.root", "RECREATE");

	Pythia8::Pythia pythia;
	pythia.readString("HardQCD:all = on");
	pythia.readString("PhaseSpace:pTHatMin = 20.");
	pythia.readString("Beams:eCM = 200.");

	pythia.init();

	std::vector<fourVec> gammas;

	std::vector<fourVec> etaGammas1;
	std::vector<fourVec> etaGammas2;

	std::vector<fourVec> piPluses;
	std::vector<fourVec> piMinuses;
	//std::vector<fourVec> piZeroes;

	fourVec temp;

	TRandom3 *rand = new TRandom3();

	TH1D *gammasHist = new TH1D("gammasHist", "gammasHist", 500,0,900);
	TH1D *etaPrimeHist = new TH1D("etaPrimeHist", "etaPrimeHist", 700,800,1100);
	TH1D *etaPrimeHistNew = new TH1D("etaPrimeHistNew", "etaPrimeHistNew", 700,900,1000);

	for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

  		if (!pythia.next()) continue;

    		for (int i = 0; i < pythia.event.size(); ++i) {
			
			if (pythia.event[i].isFinal() && pythia.event[i].id() == 22) {

				temp.E  = rand->Gaus(pythia.event[i].e(),  gSMEARING_FACTOR);		
				temp.px = rand->Gaus(pythia.event[i].px(), gSMEARING_FACTOR);		
				temp.py = rand->Gaus(pythia.event[i].py(), gSMEARING_FACTOR);		
				temp.pz = rand->Gaus(pythia.event[i].pz(), gSMEARING_FACTOR);		

				gammas.push_back( temp );
			}

			if (pythia.event[i].isFinal() && pythia.event[i].id() == 211) {

				temp.E  = rand->Gaus(pythia.event[i].e(),  gSMEARING_FACTOR);		
				temp.px = rand->Gaus(pythia.event[i].px(), gSMEARING_FACTOR);		
				temp.py = rand->Gaus(pythia.event[i].py(), gSMEARING_FACTOR);		
				temp.pz = rand->Gaus(pythia.event[i].pz(), gSMEARING_FACTOR);		

				piPluses.push_back( temp );
			}

			if (pythia.event[i].isFinal() && pythia.event[i].id() == -211) {

				temp.E  = rand->Gaus(pythia.event[i].e(),  gSMEARING_FACTOR);		
				temp.px = rand->Gaus(pythia.event[i].px(), gSMEARING_FACTOR);		
				temp.py = rand->Gaus(pythia.event[i].py(), gSMEARING_FACTOR);		
				temp.pz = rand->Gaus(pythia.event[i].pz(), gSMEARING_FACTOR);		

				piMinuses.push_back( temp );
			}
			
			/*if (pythia.event[i].id() == 221) {

				temp.E  = rand->Gaus(pythia.event[i].e(),  gSMEARING_FACTOR);		
				temp.px = rand->Gaus(pythia.event[i].px(), gSMEARING_FACTOR);		
				temp.py = rand->Gaus(pythia.event[i].py(), gSMEARING_FACTOR);		
				temp.pz = rand->Gaus(pythia.event[i].pz(), gSMEARING_FACTOR);		

				piZeroes.push_back( temp );
			}*/

		}
		
		FindEtaGammas( gammas, etaGammas1, etaGammas2, gammasHist );

		MixParticle( etaGammas1, etaGammas2, piPluses, piMinuses, etaPrimeHist );

		//MixParticleNew( piZeroes, piPluses, piMinuses, etaPrimeHistNew );

		gammas.clear();

		etaGammas1.clear();
		etaGammas2.clear();

		piPluses.clear();
		piMinuses.clear();
		//piZeroes.clear();

	}

	gammasHist->Write();
	etaPrimeHist->Write();
	etaPrimeHistNew->Write();
	
	outfile->Close();

 	pythia.stat();

	return 0;
}


void FindEtaGammas( std::vector<fourVec> &particles,
	            std::vector<fourVec> &etaGammas1,
	            std::vector<fourVec> &etaGammas2, 
		    TH1D *hist) {

	double M;
	double E1, px1, py1, pz1;
	double E2, px2, py2, pz2;

	for (unsigned int iParticle = 0; iParticle < particles.size(); iParticle++) {
		for (unsigned int jParticle = iParticle + 1; jParticle < particles.size(); jParticle++) {

			E1  = particles[iParticle].E;
			px1 = particles[iParticle].px;
			py1 = particles[iParticle].py;
			pz1 = particles[iParticle].pz;

			E2  = particles[jParticle].E;
			px2 = particles[jParticle].px;
			py2 = particles[jParticle].py;
			pz2 = particles[jParticle].pz;

			M = 1000*sqrt( (E1+E2)   * (E1+E2) 
				 -(px1+px2) * (px1+px2)
				 -(py1+py2) * (py1+py2)
				 -(pz1+pz2) * (pz1+pz2) );

			hist->Fill(M);

			if ( fabs(M - gETA_MASS) < 3 * gETA_SIGMA ) {
				etaGammas1.push_back( particles[iParticle] );	
				etaGammas2.push_back( particles[jParticle] );	
			}
		}
	}	

	return;
}

void MixParticle( std::vector<fourVec> &etaGammas1, 
		  std::vector<fourVec> &etaGammas2, 
		  std::vector<fourVec> &piPluses, 
		  std::vector<fourVec> &piMinuses, 
		  TH1D *hist ) {
	
	double M;
	double E1, px1, py1, pz1;
	double E2, px2, py2, pz2;
	double E3, px3, py3, pz3;
	double E4, px4, py4, pz4;
	
	for (unsigned int iGamma = 0; iGamma < etaGammas2.size(); iGamma++) {
		for (unsigned int iPi = 0; iPi < piPluses.size(); iPi++) {
			for (unsigned int jPi = 0; jPi < piMinuses.size(); jPi++) {

				E1  = etaGammas1[iGamma].E;
				px1 = etaGammas1[iGamma].px;
				py1 = etaGammas1[iGamma].py;
				pz1 = etaGammas1[iGamma].pz;

				E2  = etaGammas2[iGamma].E;
				px2 = etaGammas2[iGamma].px;
				py2 = etaGammas2[iGamma].py;
				pz2 = etaGammas2[iGamma].pz;

				E3  = piPluses[iPi].E;
				px3 = piPluses[iPi].px;
				py3 = piPluses[iPi].py;
				pz3 = piPluses[iPi].pz;

				E4  = piMinuses[jPi].E;
				px4 = piMinuses[jPi].px;
				py4 = piMinuses[jPi].py;
				pz4 = piMinuses[jPi].pz;

				M = 1000*sqrt( (E1+E2+E3+E4)   * (E1+E2+E3+E4) 
					 -(px1+px2+px3+px4) * (px1+px2+px3+px4)
					 -(py1+py2+py3+py4) * (py1+py2+py3+py4)
					 -(pz1+pz2+pz3+pz4) * (pz1+pz2+pz3+pz4) );

				hist->Fill(M);

			}
		}
	}

	return;
}	

void MixParticleNew( std::vector<fourVec> &piZeroes, 
		  std::vector<fourVec> &piPluses, 
		  std::vector<fourVec> &piMinuses, 
		  TH1D *hist ) {
	
	double M;
	double E1, px1, py1, pz1;
	double E2, px2, py2, pz2;
	double E3, px3, py3, pz3;
	
	for (unsigned int iGamma = 0; iGamma < piZeroes.size(); iGamma++) {
		for (unsigned int iPi = 0; iPi < piPluses.size(); iPi++) {
			for (unsigned int jPi = 0; jPi < piMinuses.size(); jPi++) {

				E1  = piZeroes[iGamma].E;
				px1 = piZeroes[iGamma].px;
				py1 = piZeroes[iGamma].py;
				pz1 = piZeroes[iGamma].pz;

				E2  = piPluses[iPi].E;
				px2 = piPluses[iPi].px;
				py2 = piPluses[iPi].py;
				pz2 = piPluses[iPi].pz;

				E3  = piMinuses[jPi].E;
				px3 = piMinuses[jPi].px;
				py3 = piMinuses[jPi].py;
				pz3 = piMinuses[jPi].pz;

				M = 1000*sqrt( (E1+E2+E3)   * (E1+E2+E3) 
					 -(px1+px2+px3) * (px1+px2+px3)
					 -(py1+py2+py3) * (py1+py2+py3)
					 -(pz1+pz2+pz3) * (pz1+pz2+pz3) );

				hist->Fill(M);

			}
		}
	}

	return;
}	



















