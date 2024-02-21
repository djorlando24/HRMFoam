   double Eth_rhol, Gas_rhol, Original_rho_l;
   
   if (Fuel.Flag_eos() ==6.0)
   {

     Eth_rhol = Rackett(1.0, 0.0, Fuel.Tcrit_a, Fuel.Tcrit_g, Fuel.Pcrit_a, Fuel.Pcrit_g, Fuel.Vcrit_a, Fuel.Vcrit_g,
                     Fuel.ZRA_a, Fuel.ZRA_g, Tref);      
     Gas_rhol = Rackett(0.0, 1.0, Fuel.Tcrit_a, Fuel.Tcrit_g, Fuel.Pcrit_a, Fuel.Pcrit_g, Fuel.Vcrit_a, Fuel.Vcrit_g,
                     Fuel.ZRA_a, Fuel.ZRA_g, Tref);     
     Original_rho_l = Rackett(Mol_alc, Mol_gas, Fuel.Tcrit_a, Fuel.Tcrit_g, Fuel.Pcrit_a, Fuel.Pcrit_g, Fuel.Vcrit_a, Fuel.Vcrit_g,
                     Fuel.ZRA_a, Fuel.ZRA_g, Tref); 
 }
   else
   {  
     EqOfSt Ethanol_1(Fuel.Tcrit_a, Fuel.Pcrit_a, Tref, Fuel.ActFt_a, Fuel.Flag_eos(),  Fuel.K1_a);
     EqOfSt Gasoline_1(Fuel.Tcrit_g, Fuel.Pcrit_g, Tref, Fuel.ActFt_g, Fuel.Flag_eos(), Fuel.K1_g);
     MixRule Ethanol_2(Gasahol_orig.P_alcohol, Tref, 1.0, 0.0, Ethanol_1.a_eos(), Ethanol_1.b_eos(), Gasoline_1.a_eos(), Gasoline_1.b_eos(),
                     Ethanol_1.u_eos(), Ethanol_1.w_eos());
     Eth_rhol = Ethanol_2.SatDen_liq;

     MixRule Gasoline_2(Gasahol_orig.P_gasoline, Tref, 0.0, 1.0, Ethanol_1.a_eos(), Ethanol_1.b_eos(), Gasoline_1.a_eos(), Gasoline_1.b_eos(),
                      Gasoline_1.u_eos(), Gasoline_1.w_eos());  
    Gas_rhol = Gasoline_2.SatDen_liq; 

    MixRule Gasahol_orig2(Gasahol_orig.PVap, Tref, Mol_alc, Mol_gas, Ethanol_1.a_eos(), Ethanol_1.b_eos(), Gasoline_1.a_eos(),
                        Gasoline_1.b_eos(),Ethanol_1.u_eos(), Ethanol_1.w_eos());
 
    Original_rho_l = Gasahol_orig2.SatDen_liq ; 

   }

//converting rhol to kg/m^3
    Eth_rhol = Eth_rhol* Fuel.MolWt_a *1E-3 ;
      
    Gas_rhol = Gas_rhol* Fuel.MolWt_g *1E-3 ;
   
    Original_rho_l = Original_rho_l * Fuel.MixMol*1E-3;

    double MassFract = Mol_alc* Fuel.MolWt_a/(Mol_alc* Fuel.MolWt_a + Mol_gas* Fuel.MolWt_g);

//      Eth_rhol = 790.00;
//     Gas_rhol = 767.00;
  
    double VolFract  = MassFract /Eth_rhol /(MassFract/Eth_rhol + (1-MassFract)/Gas_rhol); 

   cout<<"Ethanol Density (kg/m^3) "<<Eth_rhol   <<endl;
   cout<<"Volume Fraction of Ethanol " <<VolFract*100<<" % "<<endl;
   cout<<"Mass Fraction of Ethanol " <<MassFract*100<<" % "<<endl;

   
