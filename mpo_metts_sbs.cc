#include "itensor/all.h"
#include "itensor/util/stats.h"
//#include "hams/Heisenberg.h"

using namespace itensor;
using std::vector;
using namespace std;

//char response;
//std::cin>>response;


Real c1(int m, Real M){
  auto x = m/M;
  auto c1y = 0.5*(1-tanh((x-0.5)/(x*(1-x))));
  return c1y;
  //return 1.0;
}

vector<ITensor>
makeH0(SiteSet const& sites)
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
makeH(SiteSet const& sites)
    {
    auto ti = 0.0;
    auto M = 10.0;
    auto N = sites.N();
    auto H = vector<ITensor>(N+1);
    for(auto b : range1(N-1))
        {
        //Make S.S operator on sites b, b+1
        if (b>= 1 && b<=M) {
          ti = c1(M-b,M); }
        else if (b>M && b<=N-M){
          ti = 1;
        } else if(b>N-M && b<N){
          ti = c1(b-N+M,M);
        }

        H.at(b) = ti*sites.op("Sz",b)*sites.op("Sz",b+1)
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

    int N = 50;
    Real tau = 0.05;
    Real cutoff = 1E-12;
    int nsweep = 10,heatstep = 10;
    auto Energy = 0.0,Energysum = 0.0,ensemble_sum = 0.0,Finite_energy = 0.0,Energy2 = 0.0,variance = 0.0,Eave = 0.0;
    Real get_rand, Pszu,Pszd,Psxu,Psxd;
    Real cps[N],FiniteTE[20];
    Real Energybond[N],Energybond_1[N],qnbond[N+1],qnbondsum[N+1],Energybondsum[N],Energybondsum_1[N],rsz,rszsum;

    //Metts parametor
    int thermal_step = 1049,ii = 0, snum = 10;
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

    //auto H = makeH(sites);
    auto effH = makeH(sites);
    auto effH0 = makeH0(sites);
    //auto gates = makeGates(H,tau/2.0);
    //auto gates = makeGates(H,tau);

    auto ampo = AutoMPO(sites);
    int M = 10;
    auto ti = 0.0;
    //Make the Heisenberg Hamiltonian
    for(int b = 1; b < N; ++b)
        {
          if (b>= 1 && b<=M) {
            ti = c1(M-b,M); }
          else if (b>M && b<=N-M){
            ti = 1;
          } else if(b>N-M && b<N){
            ti = c1(b-N+M,M);
          }
          ampo += ti*0.5,"S+",b,"S-",b+1;
          ampo += ti*0.5,"S-",b,"S+",b+1;
          ampo +=  ti,   "Sz",b,"Sz",b+1;
        }
    auto H = MPO(ampo);

    ofstream outfile_FiniteE,outfile_config; //outfile_C;outfile_Sn,outfile_SS;
    ofstream outfile_bond,outfile_bond_1,outfile_qnbond,outfile_qnbond_1;
    //outfile_ensemble.open("ensemble.txt"); outfile_C.open("myfile_C.txt"); outfile_Sn.open("myfile_Sn.txt"); outfile_SS.open("myfile_SS.txt");
    outfile_FiniteE.open("Metts_FE.txt");
    outfile_config.open("configuration.txt");
    outfile_bond_1.open("general_bond_1.txt");
    outfile_bond.open("general_bond.txt");
    outfile_qnbond.open("general_qn.txt");


    //PrintData(gates.at(1));

    clock_t start, end;
    srand( time(NULL) );
    auto expH = toExpH<ITensor>(ampo,tau);
    auto args = Args("Cutoff=",1E-9,"Maxm=",3000);

    for(int nt = 20; nt<=80; nt=nt+10)          {   // for different temperature

    ii = 0;
    Energysum = 0.0;
    for(int i = 0;i<=snum-1;i++) {
      binning_ave[i] = 0.0;
      binning_sum[i] = 0.0;
    }

    for(int i = 1;i < N; i++){
      Energybondsum[i] = 0.0; //measurement of energy <S.S>
    }
    for(int i = 1;i < N; i++){
      Energybondsum_1[i] = 0.0; //measurement of energy <S.S>
    }
    for(int i = 1;i <= N; i++){
      qnbondsum[i] = 0.0; //measurement of energy <S.S>
    }
    rszsum = 0.0;

  for(int i = 1; i <= N; ++i)
      {
      //Neel state
      state.set(i,i%2==1 ? "Up" : "Dn");
      }
      auto psi = MPS(state);
      for(int i = 1;i<=N;i++) {
        cps[i] = 3;
      }


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
          for(int n = 2; n <= nt; ++n)
              {
              fitApplyMPO(psi,expH,psi,args);
              //printfln("step:%2.2d Energy = %.20f",n,psiHphi(psi,effH,psi));
              }
          Energy = psiHphi(psi,H,psi);
          //printfln(" Sweep: %.2d Half 1: Energy = %.10f", tstep,Energy);

          for(int i = 1;i < N; i++){
            psi.position(i);
            auto phi = psi.A(i)*psi.A(i+1);
            auto bra = prime(phi,Site);
            Energybond_1[i] = (bra*effH.at(i)*phi).real(); //measurement of energy <S.S>
            Energybond[i] = (bra*effH0.at(i)*phi).real(); //measurement of energy <S.S>
          }

          rsz = 0.0;
          for(int i = 1;i<=N;i++) {
            psi.position(i);
            auto phi = psi.A(i);
            auto bra = prime(phi,Site);
            qnbond[i] = pow((bra*sites.op("Sz",i)*phi).real(),2); //measurement of energy <S.S>
            rsz += (bra*sites.op("Sz",i)*phi).real();
          }

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

          //outfile_config<<" "<<Energy<<" "<<tstep<<endl;
          if (tstep>=heatstep)
          {
            //Energysum += Energy;
            Energysum += Energy;
            ii += 1;
            int binnum = int(((tstep+50)/100)-1);
            binning_sum[binnum]+= Energy;
            for(int i = 1;i < N; i++){
              Energybondsum[i] += Energybond[i];
              Energybondsum_1[i] += Energybond_1[i];
            }
            for(int i = 1;i <= N; i++){
              qnbondsum[i] += qnbond[i];
            }
            rszsum += rsz;
          }

      }                   //the loop of thermal evolution


        for(int i = 1;i < N; i++){
          Energybondsum[i] = Energybondsum[i]/(thermal_step-heatstep+1);
        }

        for(int i = 1;i < N; i++){
          Energybondsum_1[i] = Energybondsum_1[i]/(thermal_step-heatstep+1);
        }

        for(int i = 1;i <= N; i++){
          qnbondsum[i] = qnbondsum[i]/(thermal_step-heatstep+1);
        }

        rsz = rszsum/(thermal_step-heatstep+1);


        for (int n =0; n<= snum-1; n++) {
          cout<<"binning energy: "<<binning_sum[n]<<endl;
          binning_ave[n] = binning_sum[n]/(100);
          cout<<n<<" binning_ave: "<<binning_ave[n]<<endl;
        }

        Finite_energy = Energysum/(thermal_step - heatstep);
        cout<<"Energysum "<<ii<<" "<<Finite_energy<<endl;
        cout<<"Temperature: "<<nt*0.1<< " E: "<<Finite_energy <<endl;

        auto binning = 0.0;
        for (int n1 = 0;n1<=snum-1;n1++) {
          binning += pow((binning_ave[n1]-Finite_energy),2);
          cout<<"n1: "<<n1<<" "<<binning<<" "<<pow((binning_ave[n1]-Finite_energy),2)<<endl;
        }
        variance = sqrt(binning/((snum-1)*snum));

        cout<<"binning_variace: "<<binning<<" "<<variance<<endl;
        outfile_FiniteE<<nt*0.1<<" "<<Finite_energy<<" "<<variance<<endl;

            for(int i = 1;i < N; i++){
                outfile_bond<<Energybondsum[i]<<" ";
                }
                outfile_bond<<endl;

            for(int i = 1;i <= N; i++){
                outfile_qnbond<<qnbondsum[i]<<" ";
                }
                outfile_qnbond<<rsz<<endl;

            for(int i = 1;i < N; i++){
                outfile_bond_1<<Energybondsum_1[i]<<" ";
                }
                outfile_bond_1<<endl;

          }

        outfile_FiniteE.close();
        outfile_config.close();
        outfile_bond.close();
        outfile_qnbond.close();
        outfile_bond_1.close();

      return 0;
      }
