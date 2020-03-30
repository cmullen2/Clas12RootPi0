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




void TestInclusive(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
 // TFile *outputFile = new TFile("/project/Gruppo3/fiber7/cmullen/outputpass1Initial/TestInvMassFile005124Pass1.root","recreate");
//  TFile *outputFile = new TFile("/project/Gruppo3/fiber7/cmullen/outputpass1Initial/V2TestInvMassFile5038Pass1.root","recreate");
  TFile *outputFile = new TFile("/project/Gruppo3/fiber7/cmullen/outputpass0/V3TestInvMassFile005124Pass0WithGammaInfo.root","recreate");
  TTree *treeVars = new TTree("ParticleVars","SomeTree");
  Double_t InvE;
  Double_t InvX;
  Double_t InvY;
  Double_t InvZ;
  Double_t InvMass;
  Double_t Mass2;
  Double_t MyInvMass;
  Double_t g1E ;
  Double_t g1X ;
  Double_t g1Y ;
  Double_t g1Z ;
  Double_t g2E ;
  Double_t g2X ;
  Double_t g2Y ;
  Double_t g2Z ;
  treeVars->Branch("InvE",&InvE);
  treeVars->Branch("InvX",&InvX);
  treeVars->Branch("InvY",&InvY);
  treeVars->Branch("InvZ",&InvZ);
  treeVars->Branch("InvMass",&InvMass);
  treeVars->Branch("Mass2",&Mass2);
  treeVars->Branch("MyInvMass",&MyInvMass);
  treeVars->Branch("g1E",&g1E);
  treeVars->Branch("g1X",&g1X);
  treeVars->Branch("g1Y",&g1Y);
  treeVars->Branch("g1Z",&g1Z);
  treeVars->Branch("g2E",&g2E);
  treeVars->Branch("g2X",&g2X);
  treeVars->Branch("g2Y",&g2Y);
  treeVars->Branch("g2Z",&g2Z);

  TChain chain;
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass1Initial/skim2_*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/pass1Initial/skim2_005124*");
  chain.Add("/project/Gruppo3/fiber7/cmullen/pass0/skim2_5124*");
//  chain.Add("/project/Gruppo3/fiber7/cmullen/RafaPass1File5038/*");
  auto files=chain.GetListOfFiles();

  //some particles
  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.6,10.6);
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
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
      
    while(c12.next()==true){

      // get particles by type
      auto electrons=c12.getByID(11); //vector of region_part_ptr
      auto gammas=c12.getByID(22);
       
      // if(electrons.size()==1 && gammas.size()==2 && protons.size()==1){
     // if(electrons.size()>0 && gammas.size()>1){
      if(electrons.size()>0 && gammas.size()>2){
       
	// set the particle momentum
	SetLorentzVector(el,electrons[0]);
	SetLorentzVector(g1,gammas[1]);
	SetLorentzVector(g2,gammas[2]);
	


	TLorentzVector inv = g1 + g2;
	InvE = inv.T();
	InvZ = inv.X();
	InvZ = inv.Y();
	InvZ = inv.Z();
	InvMass = inv.M();
	Mass2 = inv.M2();
	
	
	g1E = g1.T();
	g1X = g1.X();
	g1Y = g1.Y();
	g1Z = g1.Z();

	g2E = g2.T();
	g2X = g2.X();
	g2Y = g2.Y();
	g2Z = g2.Z();

	Double_t Init =  InvE*InvE - (InvZ*InvZ + InvY*InvY + InvX*InvX);


	if(Init<0){
		MyInvMass = -1*TMath::Sqrt(-1*Init);
}
	else{
		MyInvMass = TMath::Sqrt(Init);	
}
	



	treeVars->Fill();


	//could also get particle time etc. here too
	//Double_t eTime=electrons[0]->sci(FTOF1A)->getTime();
      }
    
       
      counter++;
    }
  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  treeVars->Write();
  treeVars->Reset();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

  outputFile->Close();

}
