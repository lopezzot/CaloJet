   Double_t        Energyem;
   Double_t        EnergyScin;
   Double_t        EnergyCher;
   Double_t        neutrinoleakage;
   Double_t        leakage;
   Double_t        NofCherenkovDetected;
   Double_t        EnergyTot;
   Double_t        PrimaryParticleEnergy;
   Char_t          PrimaryParticleName[13];
   vector<double>  *VectorSignalsR;
   vector<double>  *VectorSignalsL;
   vector<double>  *VectorSignalsCherR;
   vector<double>  *VectorSignalsCherL;
   vector<double>  *VectorL;
   vector<double>  *VectorR;

   // List of branches
   TBranch        *b_Energyem;   //!
   TBranch        *b_EnergyScin;   //!
   TBranch        *b_EnergyCher;   //!
   TBranch        *b_neutrinoleakage;   //!
   TBranch        *b_leakage;   //!
   TBranch        *b_NofCherenkovDetected;   //!
   TBranch        *b_EnergyTot;   //!
   TBranch        *b_PrimaryParticleEnergy;   //!
   TBranch        *b_PrimaryParticleName;   //!
   TBranch        *b_VectorSignalsR;   //!
   TBranch        *b_VectorSignalsL;   //!
   TBranch        *b_VectorSignalsCherR;   //!
   TBranch        *b_VectorSignalsCherL;   //!
   TBranch        *b_VectorL;   //!
   TBranch        *b_VectorR;   //!
