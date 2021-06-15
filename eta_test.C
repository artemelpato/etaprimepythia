#include <vector>
#include <cmath>

typedef std::vector<ROOT::Math::PxPyPzEVector> Event;

void generateParticle( ROOT::Math::PxPyPzEVector* );

void decayParticle( ROOT::Math::PxPyPzEVector*,
		    ROOT::Math::PxPyPzEVector*,
		    ROOT::Math::PxPyPzEVector* );

const double gMASS       = 547.862; //MeV
const double gGAMMA      = 5;       //MeV
const double gMEAN_KIN_E = 100;     //MeV

int eta_test() {

	int nEvents = 10000;

	std::vector<Event> events;

	TFile *outFile = new TFile( "output.root", "RECREATE");

	TH1D *pxHist = new TH1D("pxHist", "", 150, -700, 700);
	TH1D *pyHist = new TH1D("pyHist", "", 150, -700, 700);
	TH1D *pzHist = new TH1D("pzHist", "", 150, -700, 700);
	TH1D *enHist = new TH1D("enHist", "", 150, 500, 900);

	TH1D *gammaPxHist = new TH1D("gammaPxHist", "", 150, -700, 700);
	TH1D *gammaPyHist = new TH1D("gammaPyHist", "", 150, -700, 700);
	TH1D *gammaPzHist = new TH1D("gammaPzHist", "", 150, -700, 700);

	TH1D *invMassHist = new TH1D("invMassHist", "", 150, 400, 700);
	TH1D *bgHist =      new TH1D("bgHist",      "", 150, 400, 700);

	TH1D *multHist = new TH1D("multHist", "", 150, 0, 20);

	TRandom3 *gaus = new TRandom3();

	for(int i = 0; i < nEvents; i++) {

		int nMesons =  floor( gaus->Gaus(10, 5) );
		Event currentEvent;

		for(int j = 0; j < nMesons; j++) {

			ROOT::Math::PxPyPzEVector etaVec;
			ROOT::Math::PxPyPzEVector gammaVec1, gammaVec2;
			
			generateParticle( &etaVec );

			decayParticle( &etaVec,
				       &gammaVec1,
				       &gammaVec2 );

			currentEvent.push_back( gammaVec1 );
			currentEvent.push_back( gammaVec2 );

			pxHist->Fill( etaVec.Px() );
			pyHist->Fill( etaVec.Py() );
			pzHist->Fill( etaVec.Pz() );
			enHist->Fill( etaVec.E () );

			gammaPxHist->Fill( gammaVec1.Px() );
			gammaPyHist->Fill( gammaVec1.Py() );
			gammaPzHist->Fill( gammaVec1.Pz() );

			gammaPxHist->Fill( gammaVec2.Px() );
			gammaPyHist->Fill( gammaVec2.Py() );
			gammaPzHist->Fill( gammaVec2.Pz() );

			multHist->Fill( nMesons );
		}

		events.push_back( currentEvent );
		currentEvent.clear();
	}

	for(int i = 0; i < events.size(); i++) {
		for(int j = 0; j < events[i].size(); j++) {
			for(int k = j + 1; k < events[i].size(); k++) {

				ROOT::Math::PxPyPzEVector temp;			
				temp = events[i][j] + events[i][k];
				invMassHist->Fill( temp.M() );

			}
		}
	}

	for(int i1 = 0; i1 < events.size(); i1 += 30) {
		for(int i2 = 0; i2 < events.size(); i2 += 30) {
			for(int j = 0; j < events[i1].size(); j++) {
				for(int k = 0; k < events[i2].size(); k++) {

					ROOT::Math::PxPyPzEVector temp;			
					temp = events[i1][j] + events[i2][k];
					bgHist->Fill( temp.M() );
				}
			}
		}
	}

	TH1D *fgHist = (TH1D*) invMassHist->Clone();
	fgHist->SetName("fgHist");

	double scalingFactor = invMassHist->Integral( 100, 150 ) / 
			       bgHist->Integral( 100, 150 );

	bgHist->Scale(scalingFactor);

	fgHist->Add(bgHist, -1);
	
	pxHist->Write();
	pyHist->Write();
	pzHist->Write();
	enHist->Write();

	gammaPxHist->Write();
	gammaPyHist->Write();
	gammaPzHist->Write();

	multHist->Write();

	invMassHist->Write();
	bgHist->Write();
	fgHist->Write();

	outFile->Close();

	new TBrowser;
	return 0;
}

void generateParticle(ROOT::Math::PxPyPzEVector *fourVec) {

	static TRandom3 rand(29834);
	
	double px, py, pz, E; //four vector
	double M, T;	      //mass and kinetic energy

	M = rand.BreitWigner( gMASS, gGAMMA );
	T = rand.Exp( gMEAN_KIN_E );
	E = M + T;

	rand.Sphere( px, py, pz, sqrt( E*E - M*M ) );

	fourVec->SetPxPyPzE( px, py, pz, E );
}

void decayParticle( ROOT::Math::PxPyPzEVector* particleVec,
		    ROOT::Math::PxPyPzEVector* Vec1,
		    ROOT::Math::PxPyPzEVector* Vec2 ) {

	static TRandom3 rand(9823);

	ROOT::Math::Boost boost( particleVec->BoostToCM() );
	ROOT::Math::PxPyPzEVector boostedVec = boost( *particleVec );	

	double px, py, pz;

	rand.Sphere( px, py, pz, boostedVec.E()/2 );

	Vec1->SetPxPyPzE(  px,  py,  pz, boostedVec.E()/2 );
	Vec2->SetPxPyPzE( -px, -py, -pz, boostedVec.E()/2 );

	boost.Invert();

	*Vec1 = boost( *Vec1 );
	*Vec2 = boost( *Vec2 );
}










