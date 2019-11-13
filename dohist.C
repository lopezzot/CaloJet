/***********************************************************************
	This is the main program to process the CBNT Root Ntuple from
Athena with SUSYtup. See SUSYtup.h for more information.
	For version 7.0.0++ using name susy for Ntuple.
***********************************************************************/
#include <TROOT.h> 
#include "TTree.h" 
#include "TBranch.h" 
#include <TFile.h>
#include "MyTree.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <tuple>

using namespace std;

TFile* ftree;
MyTree bonsaiTree;

vector<fastjet::PseudoJet> inputparticles_tru;
vector<fastjet::PseudoJet> inputparticles_scin;
vector<fastjet::PseudoJet> inputparticles_cher;
vector<fastjet::PseudoJet> jetexc;
vector<fastjet::PseudoJet> jet_scin;
vector<fastjet::PseudoJet> jet_cher_t;
vector<fastjet::PseudoJet> jet_cher;
vector<fastjet::PseudoJet> jet_rec;
vector<fastjet::PseudoJet> jet_tru;
vector<TLorentzVector> muvec;
vector <double> emcomp;

vector<double> Calib_VectorScinR;
vector<double> Calib_VectorScinL;
vector<double> Calib_VectorCherR;
vector<double> Calib_VectorCherL;	

double GeV=1000.;
double pi=3.14159265;
double threshold=0.01;
double etalim=5; //5.0
//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) { 
  vector<double> calibscin(std::vector<double> vectorscin);
  vector<double> calibcher(std::vector<double> vectorcher);
  tuple<double, double, double> maptower(int index, string side);
  fastjet::PseudoJet matchjet(fastjet::PseudoJet jet_in, vector<fastjet::PseudoJet> testvec);
  fastjet::PseudoJet mergejet(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher); 
  cout << "-------------------------------------------" << std::endl;
  cout << " N arguments " << argc << std::endl;
  if(argc<2){
    cout << "Please give output.root input.root" << endl;
    return 0;
  }
  // Output tree file
  string histName=argv[1];
  std::cout << " Output file name: " << std::endl;
  std::cout << "      " << histName << std::endl;

  // Open input file
  std::string fn = argv[2];
  std::cout << " Read file name: " << std::endl;
  std::cout << "      " << fn << std::endl;

  std::string filetru = fn+"_truth.root";
  std::string filesim = fn+".root";
  std::cout<< "Reading file " <<filetru <<std::endl;
  TFile* f = new TFile( filetru.c_str() );
  std::cout<< "Reading file " <<filesim <<std::endl;
  TFile* f1 = new TFile( filesim.c_str() );

  int pos_st=histName.rfind("/");
  string stma=histName.substr(pos_st+1);
  string newfile="bonsai."+stma;
  cout << " bonsai file " << newfile << endl;
  ftree = new TFile(newfile.c_str(), "RECREATE");
  bonsaiTree.Init();
#include "truthdec.h"
#include "B4dec.h"
//
  TTree* tree1 = (TTree*)f->Get("truth");
  TTree* tree2 = (TTree*)f1->Get("B4");
  tree1->AddFriend(tree2);
#include "truthset.h"
#include "B4set.h"
  tree1->GetEntry(0);
//
  if (tree1== 0) return 1;
  Int_t nentries = Int_t(tree1->GetEntries());
  Int_t nbytes= 0, nb = 0;
//
  cout << " Number of events " << nentries << endl;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    nb = tree1->GetEntry(jentry);   nbytes += nb;
    bonsaiTree.Reset();
    inputparticles_tru.clear();
    muvec.clear();
    int nmuon=0;
    int nneu=0;
    double muene_sci=0.;
    double muene_che=0.;

    for(uint itru=0;itru<mcs_n;itru++){
      int partid = mcs_pdgId->at(itru);
      double parteta = mcs_eta->at(itru);
      if(abs(partid) != 13 &&  
         abs(partid) !=12  && abs(partid) != 14 && abs(partid) != 16 &&
         abs(partid) != 1000022){
        if(abs(parteta)<etalim){
          TLorentzVector trup;
          trup.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), mcs_phi->at(itru),
                            mcs_m->at(itru));
          fastjet::PseudoJet fj(trup.Px(), trup.Py(), trup.Pz(), trup.E());
          fj.set_user_index(itru);
          inputparticles_tru.push_back(fj);
        }
      }
      if(abs(partid) ==12  && abs(partid) == 14 && abs(partid) == 16)nneu++;
      if(abs(partid) == 13){
        TLorentzVector muon;
        muon.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), mcs_phi->at(itru),
        mcs_m->at(itru));
        muvec.push_back(muon);
        nmuon++;
      } 
    } // loop on truth particles    

    jetexc.clear();
    fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, 2.*pi, 1.);
    fastjet::ClusterSequence clust_seq(inputparticles_tru, jet_def); 
    jetexc = clust_seq.exclusive_jets(int(2));


//    cout << " phi " << jetexc[0].phi() << endl;
//
//  now the rec part
//
    Calib_VectorScinR.clear();
    Calib_VectorScinL.clear();
    Calib_VectorCherR.clear();
    Calib_VectorCherL.clear();
// 
    Calib_VectorScinR = calibscin(*VectorSignalsR);
    Calib_VectorScinL = calibscin(*VectorSignalsL);
    Calib_VectorCherR = calibcher(*VectorSignalsCherR);
    Calib_VectorCherL = calibcher(*VectorSignalsCherL);
    double energy=0;
    for(uint i=0; i<Calib_VectorScinR.size(); i++) {
      energy+=Calib_VectorScinR.at(i)+Calib_VectorScinL.at(i);
    }       
    if(energy>0){
      inputparticles_scin.clear();
      inputparticles_cher.clear();
// right side
      for(int towerindex=1; towerindex<=75*36; towerindex++) {
        auto thphieta=maptower(towerindex, "right");
        double theta=get<0>(thphieta);
        double phi=get<1>(thphieta);
        double eta=get<2>(thphieta);
//        cout << towerindex << " theta " << theta << " phi " << phi << " eta " << eta << endl;
        double energy_scin = Calib_VectorScinR[towerindex];
        double pt_scin = energy_scin*sin(theta*pi/180.);

        double energy_cher = Calib_VectorCherR[towerindex];
        double pt_cher = energy_cher*sin(theta*pi/180.);

        TLorentzVector towerscin;
        towerscin.SetPtEtaPhiM(pt_scin, eta, phi*pi/180., 0.);
        TLorentzVector towercher;
        towercher.SetPtEtaPhiM(pt_cher, eta, phi*pi/180., 0.);
        double deltamumin=999999.;
        for(uint i=0;i<muvec.size();i++) {
           double deltaR=abs(towerscin.DeltaR(muvec[i]));
           if(deltaR<deltamumin)deltamumin=deltaR;
        }
        if(energy_scin > threshold) {
          if(deltamumin<0.1){
            muene_sci = muene_sci+towerscin.E();
            muene_che = muene_che+towercher.E();
          }
          if(deltamumin>0.1){
            inputparticles_scin.push_back(fastjet::PseudoJet(towerscin.Px(), 
                                          towerscin.Py(), towerscin.Pz(), towerscin.E()));
            inputparticles_cher.push_back(fastjet::PseudoJet(towercher.Px(), 
                                          towercher.Py(), towercher.Pz(), towercher.E())); 
          }
        }  
      }
// left sid3
      for(int towerindex=1; towerindex<=75*36; towerindex++) {
        auto thphieta=maptower(towerindex, "left");
        double theta=get<0>(thphieta);
        double phi=get<1>(thphieta);
        double eta=get<2>(thphieta);

//        cout << towerindex << " theta " << theta << " phi " << phi << " eta " << eta << endl;
        double energy_scin = Calib_VectorScinL[towerindex];
        double pt_scin = energy_scin*sin(theta*pi/180.);

        double energy_cher = Calib_VectorCherL[towerindex];
        double pt_cher = energy_cher*sin(theta*pi/180.);

        TLorentzVector towerscin;
        towerscin.SetPtEtaPhiM(pt_scin, eta, phi*pi/180., 0.);
        TLorentzVector towercher;
        towercher.SetPtEtaPhiM(pt_cher, eta, phi*pi/180., 0.);
        double deltamumin=999999.;
        for(uint i=0;i<muvec.size();i++) {
           double deltaR=abs(towerscin.DeltaR(muvec[i]));
           if(deltaR<deltamumin)deltamumin=deltaR;
        }
        if(energy_scin > threshold) {
          if(deltamumin<0.1){
            muene_sci = muene_sci+towerscin.E();
            muene_che = muene_che+towercher.E();
          }
          if(deltamumin>0.1){
            inputparticles_scin.push_back(fastjet::PseudoJet(towerscin.Px(), 
                                          towerscin.Py(), towerscin.Pz(), towerscin.E()));
            inputparticles_cher.push_back(fastjet::PseudoJet(towercher.Px(), 
                                          towercher.Py(), towercher.Pz(), towercher.E())); 
          }
        }  
      }
//      cout << " input scin " << inputparticles_scin.size()  << endl;
//
      fastjet::JetDefinition jet_defs(fastjet::ee_genkt_algorithm, 2.*pi, 1.);
      fastjet::ClusterSequence clust_seq_scin(inputparticles_scin, jet_defs); 
      fastjet::ClusterSequence clust_seq_cher(inputparticles_cher, jet_defs);

//  clear vectors of jets
      jet_scin.clear();
      jet_cher.clear();
      jet_cher_t.clear();
      jet_rec.clear();
      jet_tru.clear();
//  create vector of jet_scin
      jet_scin.push_back(clust_seq_scin.exclusive_jets(int(2))[0]);
      jet_scin.push_back(clust_seq_scin.exclusive_jets(int(2))[1]);
//   create temp vector of jet_cher 
      jet_cher_t.push_back(clust_seq_cher.exclusive_jets(int(2))[0]);
      jet_cher_t.push_back(clust_seq_cher.exclusive_jets(int(2))[1]);
//   align jet_cher and jet_scin vector
      for(uint jn=0; jn<jet_scin.size();jn++) {
        jet_cher.push_back(matchjet(jet_scin[jn], jet_cher_t)); 
      }
//   combine aligned scin and cher into rec
      for(uint jn=0; jn<jet_scin.size();jn++) {
        jet_rec.push_back(mergejet(jet_scin[jn],jet_cher[jn]));
      }
//   align truth jet with rec jets
//
      for(uint jn=0; jn<jet_rec.size();jn++) {
        jet_tru.push_back(matchjet(jet_rec[jn], jetexc)); 
      }
//
//   calculate EM fraction for jets
//
    emcomp.clear();
    for(uint jt=0; jt<jet_tru.size(); jt++) {
       vector<fastjet::PseudoJet> constituents = jet_tru[jt].constituents();
       double eem=0;
       for (uint jc=0; jc<constituents.size(); jc++){
	 int nc=constituents[jc].user_index();
	 int partid = mcs_pdgId->at(nc);
	 if(abs(partid)==22 || abs(partid)==11) eem+=mcs_E->at(nc);
       }
       emcomp.push_back(eem);
    }   
	 
//
//    save in ntuple
//
      if(jet_rec.size()==2 && jet_tru.size()==2) {
        fastjet::PseudoJet jetrec=jet_rec[0]+jet_rec[1];
        fastjet::PseudoJet jettruth=jet_tru[0]+jet_tru[1];
//

        bonsaiTree.nmuon = nmuon;
        bonsaiTree.nneu = nneu;
        bonsaiTree.mjjr= jetrec.m();
        bonsaiTree.mjjt= jettruth.m();
        bonsaiTree.edep= EnergyTot/1000.;
        bonsaiTree.muene_sci=muene_sci;
        bonsaiTree.muene_che=muene_che;
	bonsaiTree.emcomp1=emcomp[0];
	bonsaiTree.emcomp2=emcomp[1];
//
        bonsaiTree.j1t_E=jet_tru[0].E();
        bonsaiTree.j1t_pt=jet_tru[0].pt();
        bonsaiTree.j1t_eta=jet_tru[0].eta();
        bonsaiTree.j1t_phi=jet_tru[0].phi();
        bonsaiTree.j1t_m=jet_tru[0].m();
        bonsaiTree.j1t_theta=jet_tru[0].theta();
        bonsaiTree.j2t_E=jet_tru[1].E();
        bonsaiTree.j2t_pt=jet_tru[1].pt();
        bonsaiTree.j2t_eta=jet_tru[1].eta();
        bonsaiTree.j2t_phi=jet_tru[1].phi();
        bonsaiTree.j2t_m=jet_tru[1].m();
        bonsaiTree.j2t_theta=jet_tru[1].theta();
        bonsaiTree.j1r_E=jet_rec[0].E();
        bonsaiTree.j1r_pt=jet_rec[0].pt();
        bonsaiTree.j1r_eta=jet_rec[0].eta();
        bonsaiTree.j1r_phi=jet_rec[0].phi();
        bonsaiTree.j1r_m=jet_rec[0].m();
        bonsaiTree.j1r_theta=jet_rec[0].theta();
        bonsaiTree.j2r_E=jet_rec[1].E();
        bonsaiTree.j2r_pt=jet_rec[1].pt();
        bonsaiTree.j2r_eta=jet_rec[1].eta();
        bonsaiTree.j2r_phi=jet_rec[1].phi();
        bonsaiTree.j2r_m=jet_rec[1].m();
        bonsaiTree.j2r_theta=jet_rec[1].theta();
        bonsaiTree.j1s_E=jet_scin[0].E();
        bonsaiTree.j1s_pt=jet_scin[0].pt();
        bonsaiTree.j1s_eta=jet_scin[0].eta();
        bonsaiTree.j1s_phi=jet_scin[0].phi();
        bonsaiTree.j1s_m=jet_scin[0].m();
        bonsaiTree.j1s_theta=jet_scin[0].theta();
        bonsaiTree.j2s_E=jet_scin[1].E();
        bonsaiTree.j2s_pt=jet_scin[1].pt();
        bonsaiTree.j2s_eta=jet_scin[1].eta();
        bonsaiTree.j2s_phi=jet_scin[1].phi();
        bonsaiTree.j2s_m=jet_scin[1].m();
        bonsaiTree.j2s_theta=jet_scin[1].theta();
        bonsaiTree.j1c_E=jet_cher[0].E();
        bonsaiTree.j1c_pt=jet_cher[0].pt();
        bonsaiTree.j1c_eta=jet_cher[0].eta();
        bonsaiTree.j1c_phi=jet_cher[0].phi();
        bonsaiTree.j1c_m=jet_cher[0].m();
        bonsaiTree.j1c_theta=jet_cher[0].theta();
        bonsaiTree.j2c_E=jet_cher[1].E();
        bonsaiTree.j2c_pt=jet_cher[1].pt();
        bonsaiTree.j2c_eta=jet_cher[1].eta();
        bonsaiTree.j2c_phi=jet_cher[1].phi();
        bonsaiTree.j2c_m=jet_cher[1].m();
        bonsaiTree.j2c_theta=jet_cher[1].theta();
      }
          
    }// energy>0          
    
// fill output tree
    bonsaiTree.Fill();
  } // loop on events
  
  ftree->cd();
  bonsaiTree.Write();
  delete f;
  delete f1;

}
std::vector<double> calibscin(std::vector<double> vectorscin){
	std::vector<double> s_cont;
	s_cont = {391.9964784225476, 392.65450934717455, 391.9130470386118, 390.6050630558208, 389.1126112802391, 387.95945028322916, 387.4823598282578, 388.1313680310084, 389.1887061971629, 389.39650703390834, 388.1539715075822, 388.50402528038387, 388.8469797524201, 389.5104813700643, 389.9263913577828, 389.03950145952115, 388.99508641977286, 388.772944352741, 389.1309191354679, 389.2350876534651, 389.13791528990845, 388.8450419930907, 389.1250633143002, 388.9702856267768, 388.89437597892635, 389.54155385650677, 389.6077230144546, 389.83471955121126, 388.9057586548166, 388.64729869880443, 389.63139958072384, 390.1488269998545, 388.91138544610874, 389.57905936367644, 389.4504769524226, 389.22625067001826, 389.79318503855836, 389.6127033215785, 389.1037017285112, 389.73580729121386, 389.31584083722316, 389.40881680859326, 389.323138549745, 389.2451419219243, 389.3922639055792, 388.796221854985, 388.97157341770327, 388.99869472269773, 389.07910256875385, 389.31098859776864, 388.7419881503485, 388.3143923895009, 389.0092636036939, 388.229773876986, 388.3854174264933, 388.05911875263905, 387.2643574001875, 387.3554637363693, 387.39402352923776, 387.02865949856425, 386.9609406993706, 386.7226638865326, 386.7789454376136, 386.7732373979207, 385.4992193760886, 385.5215522567318, 384.82993023752476, 384.24706703575663, 384.2926857786712, 383.6152585260428, 382.8533589897697, 384.67220539781823, 383.42475883749, 382.5620400263505, 376.3078422852331};
	
	int loop1 = s_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			s_cont.push_back(s_cont[i]);
		}
	}

	double c = 0.1;
	s_cont.insert(s_cont.begin(),c);

	if(s_cont.size() != vectorscin.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorScin;

	for(uint i=0; i<vectorscin.size(); i++){
		Calib_vectorScin.push_back(vectorscin[i]*(1.0/s_cont[i]));
	}

	return Calib_vectorScin;
}

std::vector<double> calibcher(std::vector<double> vectorcher){
	std::vector<double> c_cont;
	c_cont = {99.869475119436, 99.5892930698953, 99.36819945465197, 99.32851833492266, 99.23212488346117, 99.21702593539936, 99.03028402205433, 99.15912908522857, 99.2060316999247, 99.30522580785956, 99.11344505994408, 99.21283817184117, 99.25155522412936, 99.29526075724375, 99.32726514547564, 99.21535433745935, 99.24785175107958, 99.2500389422396, 99.2947777840262, 99.20923521702923, 99.30797870907622, 99.31656510767378, 99.21947678681387, 99.32713686590662, 99.32851884709454, 99.305911407397, 99.30508949248147, 99.34913871640343, 99.19261358390726, 99.16893956672706, 99.26948077177317, 99.31331547462389, 99.30086159029355, 99.37873056276479, 99.35080907470335, 99.42847084872328, 99.40510235111934, 99.30670688191663, 99.38448289299362, 99.43075203396795, 99.34555944179634, 99.38032213687576, 99.3653611014606, 99.33974923759105, 99.36938995243544, 99.41133031804752, 99.33192563010871, 99.18049429253469, 99.39705625808192, 99.24849638998961, 99.2286007870197, 99.09928229987484, 99.0993256503711, 99.19364178823666, 99.0555902591915, 99.0621809829577, 98.97675936409205, 99.13401462382535, 99.03520349760886, 99.01410173081639, 98.805439062685, 98.8931348280415, 98.84927421673683, 98.71812808136468, 98.76605765544863, 98.8150716743538, 98.664618195543, 98.52979814367153, 98.39802169218106, 98.48346797819356, 98.46914717338315, 98.21105203135059, 98.38198816987828, 98.04422762750505, 96.54297571859878};

	int loop1 = c_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			c_cont.push_back(c_cont[i]);
		}
	}

	double c = 0.1;
	c_cont.insert(c_cont.begin(),c);

	if(c_cont.size() != vectorcher.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorCher;

	for(uint i=0; i<vectorcher.size(); i++){
		Calib_vectorCher.push_back(vectorcher[i]*(1.0/c_cont[i]));
	}

	return Calib_vectorCher;
}
//
std::tuple<double, double, double> maptower(int index, string side){
//Function to return tower angles (theta and phi) given index
  int NbOfBarrel=40;
  int NbOfEndcap=35;
  int NZrot=36;
  int TotTower=NbOfBarrel+NbOfEndcap;
  index = index-1;
  int sliceindex = index/TotTower;
  int towerindex = index-(sliceindex*TotTower);
  double deltatheta = 45./(NbOfBarrel);
//get theta
  double theta = towerindex*deltatheta+deltatheta/2.;
//  cout << " thetap " << theta << endl;
//get phi
  double phi_unit = 360./NZrot;
  double phi = (sliceindex)*phi_unit;
  
  if (side == "right"){
//     cout << " thetai " << theta+90. << " phii " << phi << " etai " <<  -log(tan(((90.-theta)*pi/180./2.))) << endl;
    return std::make_tuple(theta+90., phi, -log(tan(((90.-theta)*pi/180./2.))));
  }
  if (side == "left"){
    return std::make_tuple(90.-theta, phi, log(tan(((90.-theta)*pi/180./2.))));
  }
  return std::make_tuple(0.,0.,0.);
}
fastjet::PseudoJet mergejet(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher) {

  double c=0.34;

  double jetPx = (jet_scin.px()-c*jet_cher.px())/(1-c);
  double jetPy = (jet_scin.py()-c*jet_cher.py())/(1-c);
  double jetPz = (jet_scin.pz()-c*jet_cher.pz())/(1-c);
  double jetE = (jet_scin.e()-c*jet_cher.e())/(1.-c);
  return fastjet::PseudoJet(jetPx, jetPy, jetPz, jetE);
}
fastjet::PseudoJet matchjet(fastjet::PseudoJet jet_in, vector<fastjet::PseudoJet> testvec) {

  int imin=-1;
  double deltarmin=99999.;
  for(uint i=0; i<testvec.size(); i++){
    double deltar=jet_in.delta_R(testvec.at(i));
    if(deltar<deltarmin) {
      deltarmin=deltar;
      imin=i;
    }
  }
  if(imin != -1) return testvec.at(imin);
  else
  return fastjet::PseudoJet(0., 0., 0., 0.);
} 
