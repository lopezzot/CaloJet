   // Set object pointer
   VectorSignalsR = 0;
   VectorSignalsL = 0;
   VectorSignalsCherR = 0;
   VectorSignalsCherL = 0;
   VectorL = 0;
   VectorR = 0;

   tree2->SetBranchAddress("Energyem", &Energyem, &b_Energyem);
   tree2->SetBranchAddress("EnergyScin", &EnergyScin, &b_EnergyScin);
   tree2->SetBranchAddress("EnergyCher", &EnergyCher, &b_EnergyCher);
   tree2->SetBranchAddress("NofCherenkovDetected", &NofCherenkovDetected, &b_NofCherenkovDetected);
   tree2->SetBranchAddress("neutrinoleakage", &neutrinoleakage, &b_neutrinoleakage);
   tree2->SetBranchAddress("leakage", &leakage, &b_leakage);
   tree2->SetBranchAddress("EnergyTot", &EnergyTot, &b_EnergyTot);
   tree2->SetBranchAddress("PrimaryParticleEnergy", &PrimaryParticleEnergy, &b_PrimaryParticleEnergy);
   tree2->SetBranchAddress("PrimaryParticleName", PrimaryParticleName, &b_PrimaryParticleName);
   tree2->SetBranchAddress("VectorSignalsR", &VectorSignalsR, &b_VectorSignalsR);
   tree2->SetBranchAddress("VectorSignalsL", &VectorSignalsL, &b_VectorSignalsL);
   tree2->SetBranchAddress("VectorSignalsCherR", &VectorSignalsCherR, &b_VectorSignalsCherR);
   tree2->SetBranchAddress("VectorSignalsCherL", &VectorSignalsCherL, &b_VectorSignalsCherL);
   tree2->SetBranchAddress("VectorL", &VectorL, &b_VectorL);
   tree2->SetBranchAddress("VectorR", &VectorR, &b_VectorR);
