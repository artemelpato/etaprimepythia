void draw_fit() 
{
	TFile *inFile  = new TFile( "output.root",      "READ"     );
	TFile *outFile = new TFile( "output_draw.root", "RECREATE" );
	
	TH1D *pxHist = (TH1D*) inFile->Get("pxHist");
	TH1D *pyHist = (TH1D*) inFile->Get("pyHist");
	TH1D *pzHist = (TH1D*) inFile->Get("pzHist");
	TH1D *enHist = (TH1D*) inFile->Get("enHist");

	TH1D *gammaPxHist = (TH1D*) inFile->Get("gammaPxHist");
	TH1D *gammaPyHist = (TH1D*) inFile->Get("gammaPyHist");
	TH1D *gammaPzHist = (TH1D*) inFile->Get("gammaPzHist");

	TH1D *invMassHist = (TH1D*) inFile->Get("invMassHist");
	TH1D *bgHist =      (TH1D*) inFile->Get("bgHist"     );
	TH1D *fgHist =      (TH1D*) inFile->Get("fgHist"     );

	TF1 *bw = new TF1("bw", "[0]/( (x^2 - [1]^2)^2 + [1]^2*[2]^2)", 400, 700 );

	bw->SetParameters( 12500000, 500, 5 );
	bw->SetLineWidth(2);
	bw->SetLineColor(kPink);
	bw->SetLineStyle(kDotted);

	fgHist->SetLineColor(kBlue);
	fgHist->SetLineWidth(1);

	fgHist->Fit(bw);

	invMassHist->SetLineWidth(1);
	invMassHist->SetLineColor(kBlack);

	bgHist->SetLineWidth(1);
	bgHist->SetLineColor(kRed);

	pxHist->GetYaxis()->SetTitle("Number of particles");
	pyHist->GetYaxis()->SetTitle("Number of particles");
	pzHist->GetYaxis()->SetTitle("Number of particles");
	enHist->GetYaxis()->SetTitle("Number of particles");

	pxHist->GetXaxis()->SetTitle("#it{p_{x, #eta}}, [MeV]");
	pyHist->GetXaxis()->SetTitle("#it{p_{y, #eta}}, [MeV]");
	pzHist->GetXaxis()->SetTitle("#it{p_{z, #eta}}, [MeV]");
	enHist->GetXaxis()->SetTitle("#it{E_{#eta}}, [MeV]");

	gammaPxHist->GetYaxis()->SetTitle("Number of photons");
	gammaPyHist->GetYaxis()->SetTitle("Number of photons");
	gammaPzHist->GetYaxis()->SetTitle("Number of photons");

	gammaPxHist->GetXaxis()->SetTitle("#it{p_{x, #gamma}}, [MeV]");
	gammaPyHist->GetXaxis()->SetTitle("#it{p_{y, #gamma}}, [MeV]");
	gammaPzHist->GetXaxis()->SetTitle("#it{p_{z, #gamma}}, [MeV]");

	invMassHist->GetYaxis()->SetTitle("Number of pairs");
	bgHist     ->GetYaxis()->SetTitle("Number of pairs");
	fgHist     ->GetYaxis()->SetTitle("Number of particles");

	invMassHist->GetXaxis()->SetTitle("#it{M_{#gamma#gamma}}, [MeV]");
	bgHist     ->GetXaxis()->SetTitle("#it{M_{#gamma#gamma}}, [MeV]");
	fgHist     ->GetXaxis()->SetTitle("#it{M_{#gamma#gamma}}, [MeV]");

	invMassHist->SetTitle("#it{M_{#gamma#gamma}} Distribution");
	bgHist     ->SetTitle("Backgound");
	fgHist     ->SetTitle("Foreground");

	TCanvas *etaCanvas = new TCanvas();
	etaCanvas->SetCanvasSize(4000, 1000);
	etaCanvas->Divide(4, 1);
	etaCanvas->SetName("etaCanvas");
	etaCanvas->SetTitle("#it{#eta} momenta and energy");

	etaCanvas->cd(1);
	pxHist->Draw();

	etaCanvas->cd(2);
	pyHist->Draw();

	etaCanvas->cd(3);
	pzHist->Draw();

	etaCanvas->cd(4);
	enHist->Draw();

	TCanvas *gammaCanvas = new TCanvas();
	gammaCanvas->SetCanvasSize(3000, 1000);
	gammaCanvas->Divide(3, 1);		
	gammaCanvas->SetName("gammaCanvas");
	gammaCanvas->SetTitle("#it{#gamma} momenta");

	gammaCanvas->cd(1);
	gammaPxHist->Draw();

	gammaCanvas->cd(2);
	gammaPyHist->Draw();

	gammaCanvas->cd(3);
	gammaPzHist->Draw();

	TCanvas *invCanvas = new TCanvas();
	invCanvas->SetCanvasSize(3000, 1000);
	invCanvas->Divide(3, 1);

	invCanvas->cd(1);
	invMassHist->Draw();

	invCanvas->cd(2);
	invMassHist->Draw();
	bgHist->Draw("same");

	invCanvas->cd(3);
	fgHist->Draw("e3");


	for( auto h: {pxHist, pyHist, pzHist, enHist, 
		      gammaPxHist, gammaPyHist, gammaPzHist,
		      invMassHist, bgHist, fgHist} )
       	{
		h->SetStats(kFALSE);
	}

	pxHist->Write();
	pyHist->Write();
	pzHist->Write();
	enHist->Write();

	gammaPxHist->Write();
	gammaPyHist->Write();
	gammaPzHist->Write();

	invMassHist->Write();
	bgHist     ->Write();
	fgHist     ->Write();

	etaCanvas  ->Write();
	gammaCanvas->Write();
	invCanvas  ->Write();

	//inFile->Close();
	//outFile->Close();

	new TBrowser;
}
