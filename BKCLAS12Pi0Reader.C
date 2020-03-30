#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "TDatabasePDG.h"
#include "clas12reader.h"


using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}

//Double_t DeltaTime(){return fTime - HypTime();}
//Double_t HypTime(){return fPath/HypBeta()/2.99792e+08*1E9;}
//Double_t HypBeta(){Double_t pp=fP4.P();return pp/sqrt(pp*pp+fPDGMass*fPDGMass);}

/*
  inline Float_t HS::CLAS12::CLAS12Trigger::StartTime(HS::THSParticle* part){
  //If FT shift time first
  if(part->CLAS12()->getRegion()==clas12::FT){part->ShiftTime(fTimeShiftFT);return StartTime(part->DeltaTime());}
  //else
  return StartTime(part->DeltaTime());
  }
  
  inline Float_t HS::CLAS12::CLAS12Trigger::StartTime(Float_t ptime){
  //supply chosen (e-) particle vertex time
  Float_t rftime=fEventInfo->fRFTime;
  //Find the nearest rf beam bucket
  fStartTime=fSTimePeak-4.0080160*((Int_t)(std::round(((fSTimePeak-(ptime-rftime))/4.0080160))))+rftime;
  return fStartTime;
  }
*/                


Double_t StartTime(Double_t ptime, Double_t rftime){ //Need to have RFTIME from event, particle time and two hard coded numbers 
  Double_t fBeamBucket = 4.0080160;
  Double_t fSTimePeak = 44.125;
  Double_t fStartTime = fSTimePeak - fBeamBucket*( (Int_t)(std::round( ((fSTimePeak-(ptime-rftime))/fBeamBucket)  ) )) + rftime;
  return fStartTime;
}


Double_t StartTime(clas12::region_part_ptr rp3, Double_t RF){
  if (rp3->par()->getStatus()<2000){
    //Need to shift the time by factor if it is in FT,
    Double_t FTShift = 44.125 - 44.1; //fstimepeak - ?
    return StartTime(rp3->getTime() + FTShift,RF );
  }
  else{
    return StartTime(rp3->getTime(),RF);
  }

}

/*Double_t StartTime(Double_t ptime, Double_t rftime){ //Need to have RFTIME from event, particle time and two hard cooded numbers //Why is overloading not working,declaration order
  Double_t fBeamBucket = 4.0080160;
  Double_t fSTimePeak = 44.125;
  Double_t fStartTime = fSTimePeak - fBeamBucket*( (Int_t)(std::round( ((fSTimePeak-(ptime-rftime))/fBeamBucket)  ) )) + rftime;
  return fStartTime;
  }*/


Double_t DeltaTime( clas12::region_part_ptr rp2){
  Double_t fTime= rp2->getTime();
  Double_t fPath;
  if (rp2->par()->getStatus()<2000 ){ 
    Double_t ftx = rp2->ft(FTCAL)->getX();
    Double_t fty = rp2->ft(FTCAL)->getY();
    Double_t ftz = rp2->ft(FTCAL)->getZ();
    fPath = (sqrt(ftx*ftx + fty*fty + ftz*ftz))/100;
  }
  else{
    fPath = (rp2->getPath())/100;  //Path is in cm change to m
  }
  Double_t code=rp2->getPid();
  TParticlePDG* part=TDatabasePDG::Instance()->GetParticle(code);
  Double_t fPDGMass =part->Mass() ;
  Double_t rp2P = rp2->getP();
  Double_t HypBeta= rp2P/sqrt(rp2P*rp2P + fPDGMass*fPDGMass);
  Double_t HypTime = fPath/HypBeta/2.99792e+08*1E9; //Path m, Time ns
  Double_t SomeTime = fTime - HypTime;
  return SomeTime;

}



void CLAS12Pi0Reader(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //  TString outputFile;
  TFile *outputFile = new TFile("/project/Gruppo3/fiber7/cmullen/outputpass0/Test.root","recreate");
  /////////////////////////////////////

  /////////////////////////////////////////////
  //Temporary setup of tree and branches here for testing. Automate default branch info for each particle as function at later date. First test it saves branches and can get the correct required info.

  TTree *treeVars = new TTree("ParticleVars","SomeTree");



  ////ELECTRON
  Double_t ElectronP;
  Double_t ElectronTheta;
  Double_t ElectronPhi;
  Double_t ElectronTime;
  Double_t ElectronEdep;
  Double_t ElectronDeltaE;
  Double_t ElectronRegion;
  Double_t ElectronSector;
  //Double_t ElectronStatus=0;
  Double_t ElectronDeltaTime;
  Double_t EStartTime;
  Double_t RFTime;
  Double_t ElectronShiftedTime;


  treeVars->Branch("ElectronP",&ElectronP);
  treeVars->Branch("ElectronTheta",&ElectronTheta);
  treeVars->Branch("ElectronPhi",&ElectronPhi);
  treeVars->Branch("ElectronTime",&ElectronTime);
  treeVars->Branch("ElectronEdep",&ElectronEdep);
  treeVars->Branch("ElectronDeltaE",&ElectronDeltaE);
  treeVars->Branch("ElectronRegion",&ElectronRegion);
  treeVars->Branch("ElectronSector",&ElectronSector);
  //treeVars->Branch("ElectronStatus",&ElectronStatus);
  treeVars->Branch("ElectronDeltaTime",&ElectronDeltaTime);
  treeVars->Branch("RFTime",&RFTime);
  treeVars->Branch("EStartTime",&EStartTime);
  treeVars->Branch("ElectronShiftedTime",&ElectronShiftedTime);

  ////PROTON
  Double_t ProtonP;
  Double_t ProtonTheta;
  Double_t ProtonPhi;
  Double_t ProtonTime;
  Double_t ProtonEdep;
  Double_t ProtonDeltaE;
  Double_t ProtonRegion;
  Double_t ProtonSector;
  //Double_t ProtonStatus=0;
  Double_t ProtonDeltaTime;
  Double_t ProtonShiftedTime;

  treeVars->Branch("ProtonP",&ProtonP);
  treeVars->Branch("ProtonTheta",&ProtonTheta);
  treeVars->Branch("ProtonPhi",&ProtonPhi);
  treeVars->Branch("ProtonTime",&ProtonTime);
  treeVars->Branch("ProtonEdep",&ProtonEdep);
  treeVars->Branch("ProtonDeltaE",&ProtonDeltaE);
  treeVars->Branch("ProtonRegion",&ProtonRegion);
  treeVars->Branch("ProtonSector",&ProtonSector);
  //treeVars->Branch("ProtonStatus",&ProtonStatus);
  treeVars->Branch("ProtonDeltaTime",&ProtonDeltaTime);
  treeVars->Branch("ProtonShiftedTime",&ProtonShiftedTime);

  ////GAMMA1
  Double_t Gamma1P;
  Double_t Gamma1Theta;
  Double_t Gamma1Phi;
  Double_t Gamma1Time;
  Double_t Gamma1Edep;
  Double_t Gamma1DeltaE;
  Double_t Gamma1Region;
  Double_t Gamma1Sector;
  //Double_t Gamma1Status=0;
  Double_t Gamma1DeltaTime;
  Double_t Gamma1ShiftedTime;

  treeVars->Branch("Gamma1P",&Gamma1P);
  treeVars->Branch("Gamma1Theta",&Gamma1Theta);
  treeVars->Branch("Gamma1Phi",&Gamma1Phi);
  treeVars->Branch("Gamma1Time",&Gamma1Time);
  treeVars->Branch("Gamma1Edep",&Gamma1Edep);
  treeVars->Branch("Gamma1DeltaE",&Gamma1DeltaE);
  treeVars->Branch("Gamma1Region",&Gamma1Region);
  treeVars->Branch("Gamma1Sector",&Gamma1Sector);
  //treeVars->Branch("Gamma1Status",&Gamma1Status);
  treeVars->Branch("Gamma1DeltaTime",&Gamma1DeltaTime);
  treeVars->Branch("Gamma1ShiftedTime",&Gamma1ShiftedTime);

  ////GAMMA2
  Double_t Gamma2P;
  Double_t Gamma2Theta;
  Double_t Gamma2Phi;
  Double_t Gamma2Time;
  Double_t Gamma2Edep;
  Double_t Gamma2DeltaE;
  Double_t Gamma2Region;
  Double_t Gamma2Sector;
  //Double_t Gamma2Status=0;
  Double_t Gamma2DeltaTime;
  Double_t Gamma2ShiftedTime;

  treeVars->Branch("Gamma2P",&Gamma2P);
  treeVars->Branch("Gamma2Theta",&Gamma2Theta);
  treeVars->Branch("Gamma2Phi",&Gamma2Phi);
  treeVars->Branch("Gamma2Time",&Gamma2Time);
  treeVars->Branch("Gamma2Edep",&Gamma2Edep);
  treeVars->Branch("Gamma2DeltaE",&Gamma2DeltaE);
  treeVars->Branch("Gamma2Region",&Gamma2Region);
  treeVars->Branch("Gamma2Sector",&Gamma2Sector);
  //treeVars->Branch("Gamma2Status",&Gamma2Status);
  treeVars->Branch("Gamma2DeltaTime",&Gamma2DeltaTime);
  treeVars->Branch("Gamma2ShiftedTime",&Gamma2ShiftedTime);



  ////////////////////////////////////////////

  TChain chain;
  //  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5038*");
  //get the hipo data
  //   reader.open(inputFile.Data());
  auto files=chain.GetListOfFiles();

  //some particles
  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.6,10.6);
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector g1(0,0,0,0);
  TLorentzVector g2(0,0,0,0);

  auto* hmiss=new TH1F("missM","missM",200,-2,3);
  auto* hm2g=new TH1F("m2g","m2g",200,0,1);
  auto* hm2gCut=new TH1F("m2gCut","m2g",200,0,1);
   
  gBenchmark->Start("timer");
  int counter=0;
 
   
  for(Int_t i=0;i<files->GetEntries();i++){
    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());
    //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}
      
    //Add some event Pid based selections
    //////////c12.AddAtLeastPid(211,1); //at least 1 pi+
    //c12.addExactPid(11,1);    //exactly 1 electron
    //c12.addExactPid(211,1);    //exactly 1 pi+
    //c12.addExactPid(-211,1);    //exactly 1 pi-
    //c12.addExactPid(2212,1);    //exactly 1 proton
    //c12.addExactPid(22,2);    //exactly 2 gamma
    //////c12.addZeroOfRestPid();  //nothing else
    //////c12.useFTBased(); //and use the Pids from RECFT
      
    while(c12.next()==true){
      // c12.event()->getStartTime(); //hipo4
      // c12.head()->getStartTime(); //hipo3
      //Loop over all particles to see how to access detector info.
    
      /* 	for(auto& p : c12.getDetParticles()){
      //  get predefined selected information
      p->getTime();
      p->getDetEnergy();
      p->getDeltaEnergy();

      //check trigger bits
      //	 if(c12.checkTriggerBit(25)) cout<<"MesonExTrigger"<<endl;
      //	 else cout<<"NOT"<<endl;

      // get any detector information (if exists for this particle)
      // there should be a get function for any entry in the bank
      switch(p->getRegion()) {
      case FD :
      p->cal(PCAL)->getEnergy();
      p->cal(ECIN)->getEnergy();
      p->cal(ECOUT)->getEnergy();
      p->sci(FTOF1A)->getEnergy();
      p->sci(FTOF1B)->getEnergy();
      p->sci(FTOF2)->getEnergy();
      p->trk(DC)->getSector();
      p->che(HTCC)->getNphe();
      p->che(LTCC)->getNphe();
      //trajectories
      p->traj(LTCC)->getX();
      // p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
      break;
      case FT :
      p->ft(FTCAL)->getEnergy();
      p->ft(FTHODO)->getEnergy();
      break;
      case CD:
      p->sci(CTOF)->getEnergy();
      p->sci(CND)->getEnergy();
      break;
      }
      //   covariance matrix (comment in to see!)
      // p->covmat()->print();
      p->cmat();
      }*/

      // get particles by type
      auto electrons=c12.getByID(11); //vector of region_part_ptr
      auto gammas=c12.getByID(22);
      auto protons=c12.getByID(2212);
       
      if(electrons.size()==1 && gammas.size()==2 && protons.size()==1){
       
	// set the particle momentum
	SetLorentzVector(el,electrons[0]);
	SetLorentzVector(pr,protons[0]);
	SetLorentzVector(g1,gammas[0]);
	SetLorentzVector(g2,gammas[1]);
	
	ElectronP = el.P();
	ElectronTheta = el.Theta();
	ElectronPhi = el.Phi();
	ElectronTime = electrons[0]->getTime();
	ElectronEdep = electrons[0]->getDetEnergy();
	ElectronDeltaE = electrons[0]->getDeltaEnergy();
	ElectronRegion = electrons[0]->getRegion();
	ElectronSector = electrons[0]->getSector();
	//	ElectronStatus = el.P();
	ElectronDeltaTime = DeltaTime(electrons[0]);
        RFTime = c12.event()->getRFTime() ;
	EStartTime = StartTime(electrons[0],RFTime);
        ElectronShiftedTime = ElectronTime - EStartTime;

	ProtonP = pr.P();
	ProtonTheta = pr.Theta();
	ProtonPhi = pr.Phi();
	ProtonTime = protons[0]->getTime();
	ProtonEdep = protons[0]->getDetEnergy();
	ProtonDeltaE = protons[0]->getDeltaEnergy();
	ProtonRegion = protons[0]->getRegion();
	ProtonSector = protons[0]->getSector();
	//	ProtonStatus = pr.P();
	ProtonDeltaTime = DeltaTime(protons[0]);
        ProtonShiftedTime = ProtonTime - EStartTime;


	Gamma1P = g1.P();
	Gamma1Theta = g1.Theta();
	Gamma1Phi = g1.Phi();
	Gamma1Time = gammas[0]->getTime();
	Gamma1Edep = gammas[0]->getDetEnergy();
	Gamma1DeltaE = gammas[0]->getDeltaEnergy();
	Gamma1Region = gammas[0]->getRegion();
	Gamma1Sector = gammas[0]->getSector();
	//	Gamma1Status = g1.P();
	Gamma1DeltaTime = DeltaTime(gammas[0]);
        Gamma1ShiftedTime = Gamma1Time - EStartTime;


	Gamma2P = g2.P();
	Gamma2Theta = g2.Theta();
	Gamma2Phi = g2.Phi();
	Gamma2Time = gammas[1]->getTime();
	Gamma2Edep = gammas[1]->getDetEnergy();
	Gamma2DeltaE = gammas[1]->getDeltaEnergy();
	Gamma2Region = gammas[1]->getRegion();
	Gamma2Sector = gammas[1]->getSector();
	//	Gamma2Status = g2.P();
	Gamma2DeltaTime = DeltaTime(gammas[1]);
        Gamma2ShiftedTime = Gamma2Time - EStartTime;



	treeVars->Fill();
	TLorentzVector miss=beam+target-el-pr-g1-g2;
	hmiss->Fill(miss.M2());
	TLorentzVector pi0 = g1+g2;
	hm2g->Fill(pi0.M());
	if(TMath::Abs(miss.M2())<0.5)hm2gCut->Fill(pi0.M());

	//could also get particle time etc. here too
	//Double_t eTime=electrons[0]->sci(FTOF1A)->getTime();
      }
    
       
      counter++;
    }
  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  TCanvas* can=new TCanvas();
  can->Divide(2,1);
  can->cd(1);
  hmiss->DrawCopy();
  can->cd(2);
  hm2g->DrawCopy();
  hm2gCut->SetLineColor(2);
  hm2gCut->DrawCopy("same");
  
  hmiss->Write();
  hm2g->Write();
  hm2gCut->Write();

  treeVars->Write();
  treeVars->Reset();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

  outputFile->Close();

}
