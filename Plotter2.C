{

gStyle->SetLineWidth(2);
gStyle->SetLabelSize(0.075,"x");
gStyle->SetLabelSize(0.075,"y");
gStyle->SetTitleH(0.1);
gStyle->SetTitleW(0.6);

TCanvas c1;
//c1.Divide(2,4);
c1.cd(1);
ParticleVars->Draw("MissMass2>>h1(100,-0.3,0.3)","Gamma1P>0.5 && Gamma2P>0.5 && InvMass>0.12 && InvMass<0.15  ");
h1->SetLineWidth(2);
h1->SetTitle("MissMass2");
c1.Draw();
c1.SaveAs("AExclusiveClas12root5038MissMass2.png");


TCanvas c2;
//c1.cd(2);
ParticleVars->Draw("InvMass>>h2(100,0,0.3)","Gamma1P>0.5 && Gamma2P>0.5");
h2->SetLineWidth(2);
h2->SetTitle("InvMass");
c2.Draw();
c2.SaveAs("AExclusiveClas12root5038InvMass.png");

TCanvas c3;
ParticleVars->Draw("MissMassEPi>>h3(100,-5,5)","Gamma1P>0.5 && Gamma2P>0.5 && InvMass>0.12 && InvMass<0.15");
h3->SetLineWidth(2);
h3->SetTitle("MissMassEPi");
c3.Draw();
c3.SaveAs("AExclusiveClas12root5038MissMassEPi.png");
//c3.SaveAs("FTFTCDExclusiveMissMassEPiWithCuts.gif");
//c3.SaveAs("FTFTCDExclusiveMissMassEPiWithCuts.pdf");

TCanvas c4;
//c1.cd(4);
ParticleVars->Draw("MissMassEP2>>h4(100,-5,5)","Gamma1P>0.5 && Gamma2P>0.5 && InvMass>0.12 && InvMass<0.15 ");
h4->SetLineWidth(2);
h4->SetTitle("MissMassEP2");
c4.Draw();
c4.SaveAs("AExclusiveClas12root5038MissMassEP2.png");
//c4.SaveAs("FTFTCDExclusiveMissMassEP2WithCuts.gif");
//c4.SaveAs("FTFTCDExclusiveMissMassEP2WithCuts.pdf");



}



