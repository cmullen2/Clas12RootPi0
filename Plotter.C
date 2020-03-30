{

gStyle->SetLineWidth(2);
gStyle->SetLabelSize(0.075,"x");
gStyle->SetLabelSize(0.075,"y");
gStyle->SetTitleH(0.1);
gStyle->SetTitleW(0.6);

TCanvas c1;
//c1.Divide(2,4);
c1.cd(1);
ParticleVars->Draw("MissMass2>>h1(100,-0.3,0.3)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1  && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5  && TMath::Abs(ProtonDeltaTime)<1");
//ParticleVars->Draw("MissMass2>>h1(100,-0.3,0.3)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1  && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5 && ScElE<5 && TMath::Abs(ProtonTime)<1 && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
h1->SetLineWidth(2);
h1->SetTitle("MissMass2");
c1.Draw();
c1.SaveAs("ExclusiveClas12rootAll65MissMass2.png");
//c1.SaveAs("FTFTCDExclusiveMissMass2WithCuts.gif");
//c1.SaveAs("FTFTCDExclusiveMissMass2WithCuts.pdf");


TCanvas c2;
//c1.cd(2);
ParticleVars->Draw("InvMass>>h2(100,0,0.3)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5  && t>-0.5 && TMath::Abs(ProtonDeltaTime)<1");
//ParticleVars->Draw("InvMass>>h2(100,0,0.3)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5  && t>-0.5 && ScElE<5 && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
h2->SetLineWidth(2);
h2->SetTitle("InvMass");
c2.Draw();
c2.SaveAs("ExclusiveClas12rootAll65InvMass.png");
//c2.SaveAs("FTFTCDExclusiveInvMassWithCuts.gif");
//c2.SaveAs("FTFTCDExclusiveInvMassWithCuts.pdf");

TCanvas c3;
//c1.cd(3);
//ParticleVars->Draw("MissMassEPi>>h3(100,-5,5)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && InvMass>0.12 && InvMass<0.15 && t>-0.5 && ScElE<5 && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
ParticleVars->Draw("MissMassEPi>>h3(100,-5,5)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && InvMass>0.12 && InvMass<0.15 && t>-0.5  && TMath::Abs(ProtonDeltaTime)<1  ");
h3->SetLineWidth(2);
h3->SetTitle("MissMassEPi");
c3.Draw();
c3.SaveAs("ExclusiveClas12rootAll65MissMassEPi.png");
//c3.SaveAs("FTFTCDExclusiveMissMassEPiWithCuts.gif");
//c3.SaveAs("FTFTCDExclusiveMissMassEPiWithCuts.pdf");

TCanvas c4;
//c1.cd(4);
ParticleVars->Draw("MissMassEP2>>h4(100,-5,5)","Gamma1P>0.5 && Gamma2P>0.1  &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5  && TMath::Abs(ProtonDeltaTime)<1 ");
//ParticleVars->Draw("MissMassEP2>>h4(100,-5,5)","Gamma1P>0.5 && Gamma2P>0.1  &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5 && ScElE<5  && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000 ");
h4->SetLineWidth(2);
h4->SetTitle("MissMassEP2");
c4.Draw();
c4.SaveAs("ExclusiveClas12rootAll65MissMassEP2.png");
//c4.SaveAs("FTFTCDExclusiveMissMassEP2WithCuts.gif");
//c4.SaveAs("FTFTCDExclusiveMissMassEP2WithCuts.pdf");

TCanvas c5;
//c1.cd(5);
ParticleVars->Draw("CMPhi>>h5(100,-180,180)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5 &&  TMath::Abs(ProtonDeltaTime)<1  ");
//ParticleVars->Draw("CMPhi>>h5(100,-180,180)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5 && ScElE<5 && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
h5->SetLineWidth(2);
h5->SetTitle("CMPhi");
c5.Draw();
c5.SaveAs("ExclusiveClas12rootAll65CMPhi.png");
//c5.SaveAs("FTFTCDExclusiveCMPhiWithCuts.gif");
//c5.SaveAs("FTFTCDExclusiveCMPhiWithCuts.pdf");

TCanvas c6;
//c1.cd(6);
ParticleVars->Draw("W>>h6(100,1,10)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5 &&  TMath::Abs(ProtonDeltaTime)<1");
//ParticleVars->Draw("W>>h6(100,1,10)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5 && ScElE<5 && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
h6->SetLineWidth(2);
h6->SetTitle("W");
c6.Draw();
c6.SaveAs("ExclusiveClas12rootAll65WWithCuts.png");
//c6.SaveAs("FTFTCDExclusiveWWithCuts.gif");
//c6.SaveAs("FTFTCDExclusiveWWithCuts.pdf");

TCanvas c7;
//c1.cd(7);
ParticleVars->Draw("-t>>h7(100,0,1)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15  && TMath::Abs(ProtonDeltaTime)<1");
//ParticleVars->Draw("-t>>h7(100,0,1)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15  && ScElE<5 && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
h7->SetLineWidth(2);
h7->SetTitle("-t");
c7.Draw();
c7.SaveAs("ExclusiveClas12rootAll65tWithCuts.png");
//c7.SaveAs("FTFTCDExclusivetWithCuts.gif");
//c7.SaveAs("FTFTCDExclusivetWithCuts.pdf");

/*
TCanvas c8;
//c1.cd(8);
ParticleVars->Draw("ScElE>>h8(100,0,10)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5  && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
h8->SetLineWidth(2);
h8->SetTitle("ScElE");
c8.Draw();
c8.SaveAs("FTFTCDExclusiveScElEWithCuts.gif");
c8.SaveAs("FTFTCDExclusiveScElEWithCuts.pdf");
*/
TCanvas c9;
//c1.cd(5);
ParticleVars->Draw("ProtonTheta>>h9(100,0,180)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5  && TMath::Abs(ProtonDeltaTime)<1  ");
//ParticleVars->Draw("ProtonTheta>>h9(100,0,180)","Gamma1P>0.5 && Gamma2P>0.1 && TMath::Abs(MissMassEP2)<1 &&   TMath::Abs(MissMass2)<0.1 && MissMassEPi<1.5 && MissMassEPi>0.5 && InvMass>0.12 && InvMass<0.15 && t>-0.5 && ScElE<5 && TMath::Abs(ProtonTime)<1  && Gamma1Status<2000 && Gamma2Status<2000 && ElectronStatus<2000");
h9->SetLineWidth(2);
h9->SetTitle("ProtonTheta");
c9.Draw();
c9.SaveAs("ExclusiveClas12rootAll65CMPhi.png");
//c9.SaveAs("FTFTCDExclusiveCMPhiWithCuts.gif");
//c9.SaveAs("FTFTCDExclusiveCMPhiWithCuts.pdf");



}



