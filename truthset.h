// Set object pointer
   mcs_E = 0;
   mcs_pt = 0;
   mcs_m = 0;
   mcs_eta = 0;
   mcs_phi = 0;
   mcs_status = 0;
   mcs_barcode = 0;
   mcs_pdgId = 0;
   mcs_charge = 0;
   mcs_vx_x = 0;
   mcs_vx_y = 0;
   mcs_vx_z = 0;

   tree1->SetBranchAddress("mcs_n", &mcs_n, &b_mcs_n);
   tree1->SetBranchAddress("mcs_E", &mcs_E, &b_mcs_E);
   tree1->SetBranchAddress("mcs_pt", &mcs_pt, &b_mcs_pt);
   tree1->SetBranchAddress("mcs_m", &mcs_m, &b_mcs_m);
   tree1->SetBranchAddress("mcs_eta", &mcs_eta, &b_mcs_eta);
   tree1->SetBranchAddress("mcs_phi", &mcs_phi, &b_mcs_phi);
   tree1->SetBranchAddress("mcs_status", &mcs_status, &b_mcs_status);
   tree1->SetBranchAddress("mcs_barcode", &mcs_barcode, &b_mcs_barcode);
   tree1->SetBranchAddress("mcs_pdgId", &mcs_pdgId, &b_mcs_pdgId);
   tree1->SetBranchAddress("mcs_charge", &mcs_charge, &b_mcs_charge);
   tree1->SetBranchAddress("mcs_vx_x", &mcs_vx_x, &b_mcs_vx_x);
   tree1->SetBranchAddress("mcs_vx_y", &mcs_vx_y, &b_mcs_vx_y);
   tree1->SetBranchAddress("mcs_vx_z", &mcs_vx_z, &b_mcs_vx_z);
