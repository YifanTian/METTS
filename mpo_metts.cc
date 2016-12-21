#include "itensor/all.h"
#include "itensor/util/stats.h"
//#include "hams/Heisenberg.h"

using namespace itensor;
using std::vector;
using namespace std;
//std::cin>>response;

char response;

vector<ITensor>
makeH(SiteSet const& sites)
    {
    auto N = sites.N();
    auto H = vector<ITensor>(N+1);
    for(auto b : range1(N-1))
        {
        //Make S.S operator on sites b, b+1
        H.at(b) = sites.op("Sz",b)*sites.op("Sz",b+1)
                   + 0.5*sites.op("S+",b)*sites.op("S-",b+1)
                   + 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        }
    return H;
    }

vector<ITensor>
makeGates(vector<ITensor> const& H,
          Real tau)
    {
    auto gates = H;
    for(auto& g : gates)
        {
        if(!g) continue;
        g = expHermitian(-tau*g);
        }
    return gates;
    }

void
doSVD(MPS & psi,
      int b,
      ITensor phi,
      Direction dir,
      Real cutoff)
    {
    auto U = psi.A(b);
    ITensor D,V;
    svd(phi,U,D,V,{"Cutoff",cutoff});
    if(dir == Fromleft)
        {
        //multiply D into V
        psi.setA(b,U);
        psi.setA(b+1,D*V);
        }
    else
        {
        //multiply D into U
        psi.setA(b,U*D);
        psi.setA(b+1,V);
        }
    }

ITensor
applyGate(ITensor phi,
          ITensor gate)
    {
    //TODO:

    phi = phi*gate;     //1. Apply gate to phi using * operator
    phi.mapprime(1,0);      //2. Restore original prime level of phi's indices
    phi /= norm(phi);      //3. Normalize phi by dividing by norm(phi)

    return phi;
    }


int
main()
    {
    char response;

    int N = 10;
    Real tau = 0.05;
    Real cutoff = 1E-12;
    int nsweep = 10,heatstep = 10;
    auto Energybond = 0.0;
    auto Energy = 0.0,Energysum = 0.0,ensemble_sum = 0.0,Finite_energy = 0.0,Energy2 = 0.0,variance = 0.0,Eave = 0.0;
    Real get_rand, Pszu,Pszd,Psxu,Psxd;
    Real cps[N],FiniteTE[20];

    //Metts parametor
    int thermal_step = 200,ii = 0, snum = 10;
    Real binning_ave[snum],binning_sum[snum];
    Real Etstep[thermal_step];
    for(int i = 1;i<=thermal_step;i++) Etstep[i] = 0.0;
    for(int i = 1;i<=snum;i++)
    {
      binning_ave[i] = 0.0;
      binning_sum[i] = 0.0;
    }

    auto sites = SpinHalf(N);
    auto state = InitState(sites);          //??

    auto H = makeH(sites);
    //auto gates = makeGates(H,tau/2.0);
    //auto gates = makeGates(H,tau);

    auto ampo = AutoMPO(sites);
    //Make the Heisenberg Hamiltonian
    for(int b = 1; b < N; ++b)
        {
        ampo += 0.5,"S+",b,"S-",b+1;
        ampo += 0.5,"S-",b,"S+",b+1;
        ampo +=     "Sz",b,"Sz",b+1;
        }
    auto effH = MPO(ampo);

    ofstream outfile_thermal_2,outfile_thermal_8,outfile_ensemble,outfile_FiniteE,outfile_C,outfile_Sn,outfile_SS;
    outfile_thermal_2.open("thermal_step_2.txt"); outfile_ensemble.open("ensemble.txt"); outfile_FiniteE.open("Metts_FE.txt");
    outfile_C.open("myfile_C.txt"); outfile_Sn.open("myfile_Sn.txt"); outfile_SS.open("myfile_SS.txt");
    outfile_thermal_8.open("thermal_step_8.txt");

    //PrintData(gates.at(1));

    clock_t start, end;
    srand( time(NULL) );
    auto expH = toExpH<ITensor>(ampo,tau);
    auto args = Args("Cutoff=",1E-9,"Maxm=",3000);
    for(int t = 20;t<=80;t=t+2)          {   // for different temperature
    //int t = 200;
    //for(int m = 1;m<=1;m++)          {
    Energysum = 0.0;
    Energy2 = 0.0;
    ii = 0;

    for(int n = 1;n<=snum;n++)        // ensemble
    //for(int n = 1;n<=1;n++) {
    {
    int aft = 0;
    binning_sum[n] = 0.0;
    binning_ave[n] = 0.0;
    //Energysum = 0.0;
    auto psi = MPS(state);
    for(int tstep = 1; tstep <= thermal_step; tstep++)        //thermal_step
      {
          Energy = 0.0;
          if (tstep == 1)
          {
            for(auto n : range1(N)) {
                auto s = sites(n);
                ITensor sitesN(s);
              if ((rand()%(100-0+1)/100.0)>0.5)
                //if(n%2==1)
                {
                    sitesN.set(s(1),sqrt(2)/2);
                    sitesN.set(s(2),sqrt(2)/2);
                    psi.setA(n,sitesN);
                  }
              else  {
                    sitesN.set(s(1),sqrt(2)/2);
                    sitesN.set(s(2),-sqrt(2)/2);
                    psi.setA(n,sitesN);
                  }
              }
            }
          else {
            if(tstep%2 == 0){
                for(auto n : range1(N)) {
                  if(cps[n] == 1) state.set(n,"Up");
                  else  state.set(n,"Dn");
                }
                  psi = MPS(state);
              }
            else {
              for(auto n : range1(N)) {
                auto s = sites(n);
                ITensor sitesN(s);
                if(cps[n] == 1)
                  {
                    sitesN.set(s(1),sqrt(2)/2);
                    sitesN.set(s(2),sqrt(2)/2);
                    psi.setA(n,sitesN);
                  }
              else  {
                    sitesN.set(s(1),sqrt(2)/2);
                    sitesN.set(s(2),-sqrt(2)/2);
                    psi.setA(n,sitesN);
                  }
                }
              }
          }

          exactApplyMPO(psi,expH,psi);
          for(int n = 2; n <= t; ++n)
              {
              fitApplyMPO(psi,expH,psi,args);
              //printfln("step:%2.2d Energy = %.20f",n,psiHphi(psi,effH,psi));
              }
          Energy = psiHphi(psi,effH,psi);
          //printfln(" Sweep: %.2d Half 1: Energy = %.10f", tstep,Energy);

          //collapse
          for (auto n:range1(N))
          {
            if(n==1)  psi.position(n);
            ITensor wf = psi.Anc(n);
            auto s = sites(n);
            ITensor psim1(s),psim2(s),psimx1(s),psimx2(s),SxUp(s,prime(s)),SxDn(s,prime(s));
            psim1.set(s(1),1);psim1.set(s(2),0);
            psim2.set(s(1),0);psim2.set(s(2),1);
            psimx1.set(s(1),sqrt(2)/2);psimx1.set(s(2),sqrt(2)/2);
            psimx2.set(s(1),sqrt(2)/2);psimx2.set(s(2),-sqrt(2)/2);
            SxUp.set(s(1),prime(s)(1),0.5); SxUp.set(s(1),prime(s)(2),0.5);
            SxUp.set(s(2),prime(s)(1),0.5); SxUp.set(s(2),prime(s)(2),0.5);
            SxDn.set(s(1),prime(s)(1),0.5); SxDn.set(s(1),prime(s)(2),-0.5);
            SxDn.set(s(2),prime(s)(1),-0.5); SxDn.set(s(2),prime(s)(2),0.5);
            //PAUSE
            get_rand = (rand()%(100-0+1)/100.0);
          if(tstep%2 == 1){
            Pszu = (dag(prime(wf,Site))*sites.op("projUp",n)*wf).real();
            Pszd = (dag(prime(wf,Site))*sites.op("projDn",n)*wf).real();
            //printfln(" zup: %.10f zdn = %.10f sum = %.10f", Pszu, Pszd, Pszu+Pszd);
            if (Pszu >= get_rand)
              { cps[n]=1;
                if (n<N) psi.Anc(n+1) = (psim1*psi.Anc(n)*psi.Anc(n+1))/sqrt(Pszu); }
              else
              { cps[n]=0;
                if (n<N) psi.Anc(n+1) = (psim2*psi.Anc(n)*psi.Anc(n+1))/sqrt(Pszd); }
            }
            else {
              Psxu = (dag(prime(wf,Site))*SxUp*wf).real();
              Psxd = (dag(prime(wf,Site))*SxDn*wf).real();
              //printfln(" xup: %.10f xdn = %.10f sum = %.10f", Psxu, Psxd, Psxu+Psxd);
              if (Psxu >= get_rand)
                { cps[n]=1;
                  if (n<N) psi.Anc(n+1) = (psimx1*psi.Anc(n)*psi.Anc(n+1))/sqrt(Psxu); }
                else
                { cps[n]=0;
                  if (n<N) psi.Anc(n+1) = (psimx2*psi.Anc(n)*psi.Anc(n+1))/sqrt(Psxd); }
            }
          }

          Etstep[tstep]+=Energy;
          //cout<<tstep<<" "<<Etstep[tstep]<<" "<<Energy<<endl;

              if (tstep>heatstep)
              {
               binning_sum[n]+= Energy;
               Energysum += Energy;
               ii++;
               aft++;
               Energy2 += Energy*Energy;
               /*
               if (t == 20 && n == 1) {
                 outfile_thermal_2<<tstep<<" "<<Energy<<" "<<Energysum/ii<<" "<<sqrt(((Energy2/ii)-pow((Energysum/ii),2))/(ii-1))<<endl;
               } else if(t == 80 && n == 1) {
                 outfile_thermal_8<<tstep<<" "<<Energy<<" "<<Energysum/ii<<" "<<sqrt(((Energy2/ii)-pow((Energysum/ii),2))/(ii-1))<<endl;
               }
               */

              }

          }

            //Eave = Energysum/ii;
            //for (int step = heatstep;step<=thermal_step;step++) {
              binning_ave[n] = binning_sum[n]/aft;
              cout<<n<<" binning_ave "<<binning_ave[n]<<endl;
            //}

          }
          
          cout<< "Energysum "<<ii<<" "<<Energysum<<endl;
          //cout<< "Energy_ave "<<ii<<" "<<Energysum/ii<<endl;
          Finite_energy = Energysum/((thermal_step-heatstep)*snum);

          auto binning = 0.0;
          for (int n1 = 1;n1<=snum;n1++) {
            binning += pow((binning_ave[n1]-Finite_energy),2);
          }
          //variance = sqrt(binning/(snum*(snum-1)));
          variance = sqrt(binning/(snum-1));

          //variance = sqrt(((Energy2/ii)-(Finite_energy*Finite_energy))/(ii-1));
          //cout<<"energy2_average"<<" "<<Energy2/ii<<" "<<Finite_energy*Finite_energy<<endl;
          cout<<"binning_variace: "<<variance<<endl;

          cout<<" temperature: "<<t*0.1<< " E: "<<Finite_energy <<endl;

          outfile_FiniteE<<t*0.1<<" "<<Finite_energy<<" "<<variance<<endl;

          /*
          ensemble_sum  = 0.0;
          for(int i = heatstep+1;i<=thermal_step;i++)
          {
            Etstep[i] = Etstep[i]/10;
            ensemble_sum += Etstep[i];
            //cout<<i<<" "<<Etstep[i]<<endl;
            outfile_ensemble<<i<<" "<<Etstep[i]<<endl;
          }
          cout<< "ensemble_sum "<<ensemble_sum<<endl;
          cout<<" temperature: "<<t*0.1<< " E: "<<ensemble_sum/((thermal_step-heatstep))<<endl;
          outfile_FiniteE<<t*0.1<<" "<<ensemble_sum/((thermal_step-heatstep))<<endl;
          }
          */

        }

          outfile_thermal_2.close();
          outfile_thermal_8.close();
          outfile_ensemble.close();
          outfile_FiniteE.close();

          return 0;

    }
