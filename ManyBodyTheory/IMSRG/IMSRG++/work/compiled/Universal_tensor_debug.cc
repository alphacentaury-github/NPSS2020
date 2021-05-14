#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"

using namespace imsrg_util;

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
  Parameters PAR(argc,argv);

  string inputtbme = PAR.s("2bme");
  string input3bme = PAR.s("3bme");
  string reference = PAR.s("reference");
  string valence_space = PAR.s("valence_space");
  string basis = PAR.s("basis");
  string method = PAR.s("method");
  string flowfile = PAR.s("flowfile");
  string intfile = PAR.s("intfile");
  string core_generator = PAR.s("core_generator");
  string valence_generator = PAR.s("valence_generator");
  string fmt2 = PAR.s("fmt2");
  string denominator_delta_orbit = PAR.s("denominator_delta_orbit");
  string LECs = PAR.s("LECs");
  string scratch = PAR.s("scratch");

  int eMax = PAR.i("emax");
  int E3max = PAR.i("e3max");
  int lmax3 = PAR.i("lmax3");
  int targetMass = PAR.i("A");
  int nsteps = PAR.i("nsteps");
  int file2e1max = PAR.i("file2e1max");
  int file2e2max = PAR.i("file2e2max");
  int file2lmax = PAR.i("file2lmax");
  int file3e1max = PAR.i("file3e1max");
  int file3e2max = PAR.i("file3e2max");
  int file3e3max = PAR.i("file3e3max");

  double hw = PAR.d("hw");
  double smax = PAR.d("smax");
  double ode_tolerance = PAR.d("ode_tolerance");
  double ds_max = PAR.d("dsmax");
  double ds_0 = PAR.d("ds_0");
  double domega = PAR.d("domega");
  double omega_norm_max = PAR.d("omega_norm_max"); 
  double denominator_delta = PAR.d("denominator_delta");

  vector<string> opnames = PAR.v("Operators");

  vector<Operator> ops;



  // test 2bme file
  ifstream test(inputtbme);
  if( not test.good() )
  {
    cout << "trouble reading " << inputtbme << " exiting. " << endl;
    return 1;
  }
  test.close();
  // test 3bme file
  if (input3bme != "none")
  {
    test.open(input3bme);
    if( not test.good() )
    {
      cout << "trouble reading " << input3bme << " exiting. " << endl;
      return 1;
    }
    test.close();
  }



  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  ModelSpace modelspace = reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space);

  if (nsteps < 0)
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  modelspace.SetE3max(E3max);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);
  
  cout << "Making the operator..." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();

  Operator omegaNathan = Hbare;
  Operator RandomTensor(modelspace,2,0,0,2);
  omegaNathan.SetAntiHermitian();
  cout << "Reading in Nathan's omega" << endl;
  rw.ReadOperator_Nathan("omega_1b.dat","omega_2b.dat",omegaNathan);
//  rw.ReadOperator_Nathan("input/RANDOM_OMEGA_1b_eMax4.dat","input/RANDOM_OMEGA_2b_eMax4.dat",omegaNathan);
//  rw.ReadTensorOperator_Nathan("input/RANDOM_TENSOR_1b_eMax2.dat","input/RANDOM_TENSOR_2b_eMax2.dat",RandomTensor);
//  rw.ReadTensorOperator_Nathan("input/RANDOM_TENSOR_1b_eMax4.dat","input/RANDOM_TENSOR_2b_eMax4.dat",RandomTensor);
//  rw.ReadTensorOperator_Nathan("input/RANDOM_1b_eMax4.dat","input/RANDOM_2b_eMax4.dat",RandomTensor);
//  omegaNathan *= 0.01;

  cout << "Reading interactions..." << endl;

  #pragma omp parallel sections 
  {
    #pragma omp section
    {
    if (fmt2 == "me2j")
      rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
    else if (fmt2 == "navratil" or fmt2 == "Navratil")
      rw.ReadBareTBME_Navratil(inputtbme, Hbare);
    else if (fmt2 == "oslo" )
      rw.ReadTBME_Oslo(inputtbme, Hbare);
    else if (fmt2 == "oakridge" )
      rw.ReadTBME_OakRidge(inputtbme, Hbare);
     cout << "done reading 2N" << endl;
    }
  
    #pragma omp section
    if (Hbare.particle_rank >=3)
    {
      rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
      cout << "done reading 3N" << endl;
    }  
  }

  Hbare += Trel_Op(modelspace);

  HartreeFock hf(Hbare);
  hf.Solve();
  cout << "EHF = " << hf.EHF << endl;
  
  if (basis == "HF" and method !="HF")
    Hbare = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    Hbare = Hbare.DoNormalOrdering();

  if (method != "HF")
  {
    cout << "Perturbative estimates of gs energy:" << endl;
    double EMP2 = Hbare.GetMP2_Energy();
    cout << "EMP2 = " << EMP2 << endl; 
    double EMP3 = Hbare.GetMP3_Energy();
    cout << "EMP3 = " << EMP3 << endl; 
    cout << "To 3rd order, E = " << Hbare.ZeroBody+EMP2+EMP3 << endl;
  }

  // Calculate all the desired operators
  for (auto& opname : opnames)
  {
           if (opname == "R2_p1")        ops.emplace_back( R2_1body_Op(modelspace,"proton") );
      else if (opname == "R2_p2")        ops.emplace_back( R2_2body_Op(modelspace,"proton") );
      else if (opname == "R2_n1")        ops.emplace_back( R2_1body_Op(modelspace,"neutron") );
      else if (opname == "R2_n2")        ops.emplace_back( R2_2body_Op(modelspace,"neutron") );
      else if (opname == "Rp2")          ops.emplace_back( Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "Rn2")          ops.emplace_back( Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "Rm2")          ops.emplace_back( Rm2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "E2")           ops.emplace_back( ElectricMultipoleOp(modelspace,2) );
      else if (opname == "M1")           ops.emplace_back( MagneticMultipoleOp(modelspace,1) );
      else if (opname == "Fermi")        ops.emplace_back( AllowedFermi_Op(modelspace) );
      else if (opname == "GamowTeller")  ops.emplace_back( AllowedGamowTeller_Op(modelspace) );
      else if (opname == "R2CM")         ops.emplace_back( R2CM_Op(modelspace) );
      else if (opname == "HCM")          ops.emplace_back( HCM_Op(modelspace) );
      else if (opname == "Rso")          ops.emplace_back( RpSpinOrbitCorrection(modelspace) );
      else if (opname.substr(0,4) == "HCM_") // GetHCM with a different frequency, ie HCM_24 for hw=24
      {
         double hw_HCM;
         double hw_save = modelspace.GetHbarOmega();
         istringstream(opname.substr(4,opname.size())) >> hw_HCM;
         modelspace.SetHbarOmega(hw_HCM);
         ops.emplace_back( HCM_Op(modelspace) );
         modelspace.SetHbarOmega(hw_save);
      }
      else if (opname.substr(0,4) == "Rp2Z")
      {
        int Z_rp;
        istringstream(opname.substr(4,opname.size())) >> Z_rp;
        ops.emplace_back( Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) );
      }
      else if (opname.substr(0,4) == "Rn2Z")
      {
        int Z_rp;
        istringstream(opname.substr(4,opname.size())) >> Z_rp;
        ops.emplace_back( Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) );
      }
      else if (opname.substr(0,4) == "rhop")
      {
        double rr;
        istringstream(opname.substr(4,opname.size())) >> rr;
        ops.emplace_back( ProtonDensityAtR(modelspace,rr));
      }
      else if (opname.substr(0,4) == "rhon")
      {
        double rr;
        istringstream(opname.substr(4,opname.size())) >> rr;
        ops.emplace_back( NeutronDensityAtR(modelspace,rr));
      }
      else if (opname.substr(0,6) == "OneOcc")
      {
         map<char,int> lvals = {{'s',0},{'p',1},{'d',2},{'f',3},{'g',4},{'h',5}};
         char pn,lspec;
         int n,l,j,t;
         istringstream(opname.substr(6,1)) >> pn;
         istringstream(opname.substr(7,1)) >> n;
         istringstream(opname.substr(8,1)) >> lspec;
         istringstream(opname.substr(9,opname.size())) >> j;
         l = lvals[lspec];
         t = pn == 'p' ? -1 : 1;
         ops.emplace_back( NumberOp(modelspace,n,l,j,t) );
      }
      else if (opname.substr(0,6) == "AllOcc")
      {
         map<char,int> lvals = {{'s',0},{'p',1},{'d',2},{'f',3},{'g',4},{'h',5}};
         char pn,lspec;
         int l,j,t;
         istringstream(opname.substr(6,1)) >> pn;
         istringstream(opname.substr(7,1)) >> lspec;
         istringstream(opname.substr(8,opname.size())) >> j;
         l = lvals[lspec];
         t = pn == 'p' ? -1 : 1;
         ops.emplace_back( NumberOpAlln(modelspace,l,j,t) );
      }
      else if (opname.substr(0,9) == "protonFBC")
      {
         int nu;
         istringstream(opname.substr(9,opname.size())) >> nu;
         ops.emplace_back( FourierBesselCoeff( modelspace, nu, 8.0, modelspace.proton_orbits) );
      }
      else if (opname.substr(0,10) == "neutronFBC")
      {
         int nu;
         istringstream(opname.substr(10,opname.size())) >> nu;
         ops.emplace_back( FourierBesselCoeff( modelspace, nu, 8.0, modelspace.neutron_orbits) );
      }
      else //need to remove from the list
      {
         cout << "Unknown operator: " << opname << endl;
      }
  }


  

  for (auto& op : ops)
  {
     rw.WriteOperatorHuman(op, "bare.op");
     if (basis == "HF") op = hf.TransformToHFBasis(op);
     rw.WriteOperatorHuman(op, "HF.op");
     op = op.DoNormalOrdering();
     if (method == "MP3")
     {
       double dop = op.MP1_Eval( Hbare );
       cout << "Operator 1st order correction  " << dop << "  ->  " << op.ZeroBody + dop << endl;
     }
  }
  auto itR2p = find(opnames.begin(),opnames.end(),"Rp2");
  if (itR2p != opnames.end())
  {
    Operator& Rp2 = ops[itR2p-opnames.begin()];
    int Z = modelspace.GetTargetZ();
    int A = modelspace.GetTargetMass();
    cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
    cout << " HF charge radius = " << sqrt( Rp2.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl; 
  }
  for (int i=0;i<ops.size();++i)
  {
    Operator& op = ops[i];
    cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
  }



  
//  cout << "create Comm111 and Comm121 " << endl;
//  Operator Comm111(modelspace,2,0,0,2);
//  Operator Comm121(modelspace,2,0,0,2);
//  Operator Comm122(modelspace,2,0,0,2);

//  cout << "call comm111st..." << endl;
//  Comm111.comm111st(omegaNathan,ops[0]);
//  cout << "Write 111" << endl;
//  rw.WriteOperatorHuman(Comm111,"comm111_omega_nathan.op");
//  cout << "call comm121st..." << endl;
//  Comm121.comm121st(omegaNathan,ops[0]);
//  Comm122.comm122st(omegaNathan,ops[0]);
//  cout << "done. Now write operators" << endl;

//  rw.WriteOperatorHuman(Comm121,"comm121_omega_nathan.op");
//  rw.WriteOperatorHuman(Comm122,"comm122_omega_nathan.op");



  Operator OneComm(modelspace,2,0,0,2);
  Operator TwoComm(modelspace,2,0,0,2);
  Operator Nested122(modelspace,2,0,0,2);
  Operator Nested222pp(modelspace,2,0,0,2);
  Operator Nested121(modelspace,2,0,0,2);
  Operator Nested111(modelspace,2,0,0,2);
  Operator Nested222ph(modelspace,2,0,0,2);
  cout << "One Comm" << endl;
//  OneComm.SetToCommutator(omegaNathan,RandomTensor);
  OneComm.SetToCommutator(omegaNathan,ops[0]);
  rw.WriteOperatorHuman(OneComm,"One_comm_omega_nathan.op");
  cout << "Nested122st" << endl;
  Nested122.comm122st(omegaNathan,OneComm);
  rw.WriteOperatorHuman(Nested122,"Nested122_omega_nathan.op");
  cout << "Nested222_pp_hh_221st" << endl;
  Nested222pp.comm222_pp_hh_221st(omegaNathan,OneComm);
  rw.WriteOperatorHuman(Nested222pp,"Nested222pp_omega_nathan.op");
  cout << "Nested222ph" << endl;
  Nested222ph.comm222_phst(omegaNathan,OneComm);
  rw.WriteOperatorHuman(Nested222ph,"Nested222ph_omega_nathan.op");
  cout << "Nested121" << endl;
  Nested121.comm121st(omegaNathan,OneComm);
  rw.WriteOperatorHuman(Nested121,"Nested121_omega_nathan.op");
  cout << "Nested111" << endl;
  Nested111.comm111st(omegaNathan,OneComm);
  rw.WriteOperatorHuman(Nested111,"Nested111_omega_nathan.op");
  cout << "TwoComm" << endl;
  TwoComm.SetToCommutator(omegaNathan,OneComm);
  rw.WriteOperatorHuman(TwoComm,"Two_comm_omega_nathan.op");
//  TwoComm /= 2;
//  rw.WriteOperatorHuman(TwoComm,"Two_comm_divide2_omega_nathan.op");
  ops[0] = ops[0].BCH_Transform(omegaNathan);
  rw.WriteOperatorHuman(ops[0],"E2_full_omega_nathan.op");



/*

  Operator OneComm(modelspace,2,0,0,2);
  Operator TwoComm(modelspace,2,0,0,2);
  Operator Random122(modelspace,2,0,0,2);
  Operator Random222pp(modelspace,2,0,0,2);
  Operator Random121(modelspace,2,0,0,2);
  Operator Random111(modelspace,2,0,0,2);
  Operator Random222ph(modelspace,2,0,0,2);
  cout << "One Comm" << endl;
  OneComm.SetToCommutator(omegaNathan,RandomTensor);
  rw.WriteOperatorHuman(OneComm,"One_comm_omega_nathan_e4.op");
  cout << "Random122st" << endl;
  Random122.comm122st(omegaNathan,RandomTensor);
  rw.WriteOperatorHuman(Random122,"Random122_omega_nathan_e4.op");
  cout << "Random222_pp_hh_221st" << endl;
  Random222pp.comm222_pp_hh_221st(omegaNathan,RandomTensor);
  rw.WriteOperatorHuman(Random222pp,"Random222pp_omega_nathan_e4.op");
  cout << "Random222ph" << endl;
  Random222ph.comm222_phst(omegaNathan,RandomTensor);
  rw.WriteOperatorHuman(Random222ph,"Random222ph_omega_nathan_e4.op");
  cout << "Random121" << endl;
  Random121.comm121st(omegaNathan,RandomTensor);
  rw.WriteOperatorHuman(Random121,"Random121_omega_nathan_e4.op");
  cout << "Random111" << endl;
  Random111.comm111st(omegaNathan,RandomTensor);
  rw.WriteOperatorHuman(Random111,"Random111_omega_nathan_e4.op");
  cout << "TwoComm" << endl;
  TwoComm.SetToCommutator(omegaNathan,OneComm);
  rw.WriteOperatorHuman(TwoComm,"Two_comm_omega_nathan_e4.op");

  TwoComm = RandomTensor.BCH_Transform(omegaNathan);
  rw.WriteOperatorHuman(TwoComm,"Random_decoupled_omega_nathan_e4.op");

  Hbare = hf.GetNormalOrderedH();
  Hbare = Hbare.BCH_Transform(omegaNathan);
  rw.WriteOperatorHuman(Hbare, "H_decoupled_Nathan_e4.op");
  ops[0] = ops[0].BCH_Transform(omegaNathan);
  rw.WriteOperatorHuman(ops[0],"E2_full_omega_nathan_e4.op");
*/
 

//  CommutatorTest(omegaNathan,Hbare);
 
  if ( method == "HF" or method == "MP3")
  {
    Hbare.PrintTimes();
    return 0;
  }


  

  IMSRGSolver imsrgsolver(Hbare);
  imsrgsolver.SetReadWrite(rw);

  
  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=500;
    method = "magnus";
  }

  imsrgsolver.SetMethod(method);
  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(ds_max);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);


//   Operator Comm1(modelspace,0,0,0,2);
//   Operator Comm2(modelspace,0,0,0,2);
//   Operator trel = Trel_Op(modelspace);
//   Operator eta(imsrgsolver.GetEta() );
//   imsrgsolver.GetGenerator().Update(&Hbare,&eta);
//   CommutatorTest(eta,Hbare);
//   cout << "Now with trel..." << endl;
//   CommutatorTest(trel,Hbare);
//   Comm1.comm222_phss(eta, Hbare);
//   Comm2.comm222_pp_hh_221ss(eta, Hbare);
//   rw.WriteOperatorHuman(Comm1,"comm222phss.op");
//   rw.WriteOperatorHuman(Comm2,"comm222pp_hh_221ss.op");

  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  if (nsteps > 1) // two-step decoupling, do core first
  {
    imsrgsolver.SetGenerator(core_generator);
    imsrgsolver.Solve();
    if (method == "magnus") smax *= 2;
  }

  imsrgsolver.SetGenerator(valence_generator);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.Solve();


  // Transform all the operators
  if (method == "magnus")
  {
    Operator omega = imsrgsolver.GetOmega(0);
    rw.WriteOperatorHuman(omega,"omega.op");
    if (ops.size()>0) cout << "transforming operators" << endl;
    for (size_t i=0;i<ops.size();++i)
    {
      cout << opnames[i] << " " << flush;
      ops[i] = imsrgsolver.Transform(ops[i]);
      rw.WriteOperatorHuman(ops[i], "decoupled.op");
      cout << " (" << ops[i].ZeroBody << " ) " << endl; 
    }
    cout << endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  // If we're doing targeted normal ordering 
  // we now re-normal order wrt to the core
  // and do any remaining flow.
//  if (reference != "default"  and reference != valence_space)
  ModelSpace ms2(modelspace);
  bool renormal_order = false;
  if (modelspace.valence.size() > 0 )
  {
    renormal_order = modelspace.holes.size() != modelspace.core.size();
    if (not renormal_order)
    {
      for (auto c : modelspace.core)
      {
         if ( (find(modelspace.holes.begin(),modelspace.holes.end(),c) == modelspace.holes.end()) or (abs(1-modelspace.holes[c])>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
//  if ( modelspace.core != modelspace.holes )
  if ( renormal_order )
  {

    Hbare = imsrgsolver.GetH_s();

    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
    cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << endl;
    Hbare = Hbare.UndoNormalOrdering();

    ms2.SetReference(ms2.core); // chage the reference determinant
    Hbare.SetModelSpace(ms2);

    cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << endl;
    Hbare = Hbare.DoNormalOrdering();

    imsrgsolver.SetHin(Hbare);
    imsrgsolver.SetEtaCriterion(1e-4);
    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
    cout << "Final transformation on the operators..." << endl;
    for (auto& op : ops)
    {
      double ZeroBody_before = op.ZeroBody;
      op = op.UndoNormalOrdering();
      double ZeroBody_undo = op.ZeroBody;
      op.SetModelSpace(ms2);
      op = op.DoNormalOrdering();
      double ZeroBody_mid = op.ZeroBody;
      // transform using the remaining omegas
      op = imsrgsolver.Transform_Partial(op,nOmega);
      cout << ZeroBody_before << "   =>   " << ZeroBody_undo << "   =>   " << ZeroBody_mid<< "   =>   " << op.ZeroBody << endl;
    }
  }


  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
  if (modelspace.valence.size() > 0)
  {
    rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
    rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");

    if (method == "magnus")
    {
       for (int i=0;i<ops.size();++i)
       {
//          ops[i] = imsrgsolver.Transform(ops[i]);
          if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
          {
            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
          }
          else
          {
            rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
            rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
          }
       }
    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    cout << "Core Energy = " << setprecision(6) << imsrgsolver.GetH_s().ZeroBody << endl;
    for (int i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << endl; 
         cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl; 
      }
    }
  }




  Hbare.PrintTimes();
 
  return 0;
}

