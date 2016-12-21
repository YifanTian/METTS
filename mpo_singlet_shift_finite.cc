#include "itensor/all.h"

using namespace std;
using namespace itensor;
using std::vector;

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

int main()
{
int N = 10,upn,dnn,spinz;
auto sites = SpinHalf(N);
auto state = InitState(sites);
Real cutoff = 1E-12;
Real get_rand, Psuu,Psdd,Psud1,Psud2,Energysum = 0.0,Finite_energy = 0.0,variance = 0.0,Energy;
int thermal_step = 400,heatstep = 50,test_num = 5,ii;
Real cps[N+1];
Real binning_ave[test_num],binning_sum[test_num];

ofstream outfile_thermal_2,outfile_thermal_8,outfile_ensemble,outfile_FiniteE,outfile_config,outfile_C,outfile_Sn,outfile_SS;
outfile_thermal_2.open("thermal_step_2.txt"); outfile_ensemble.open("ensemble.txt"); outfile_FiniteE.open("Metts_FE.txt");
outfile_thermal_8.open("thermal_step_8.txt");
outfile_config.open("configuration.txt");


auto ampo = AutoMPO(sites);
//Make the Heisenberg Hamiltonian
for(int b = 1; b < N; ++b)
    {
    ampo += 0.5,"S+",b,"S-",b+1;
    ampo += 0.5,"S-",b,"S+",b+1;
    ampo +=     "Sz",b,"Sz",b+1;
    }
auto H = MPO(ampo);
auto tau = 0.05;
auto expH = toExpH<ITensor>(ampo,tau);

auto args = Args("Cutoff=",1E-9,"Maxm=",3000);
auto ttotal = 1.0;
//auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
auto nt = 10;

//for (int ent = 20; ent <= 60; ent=ent+10) {
for (int ent = 10; ent <= 80; ent=ent+2) {
  nt = ent;

  ii = 0;
  Energysum = 0.0;
  for(int i = 0;i<=test_num-1;i++) {
    binning_ave[i] = 0.0;
    binning_sum[i] = 0.0;
  }

for (int test = 1; test<=test_num; test += 1) {

for(int i = 1; i <= N; ++i)
    {
    //Neel state
    state.set(i,i%2==1 ? "Up" : "Dn");
    }
auto psi = MPS(state);
for(int i = 1;i<=N;i++) {
  cps[i] = 3;
}

//ii = 0;
//Energysum = 0.0;

//outfile_config<<"configuration: "<<"spin_z "<<" Energy"<<endl;

for(int tstep = 1; tstep <= thermal_step; tstep++)        //thermal_step
  {
    upn = 0;
    dnn = 0;
    spinz = 0;
    if (tstep%2 == 1){
      for (int c = 1; c<=N; c+=2){
        if (cps[c] == 1){
          //outfile_config<<"11 ";
          spinz += 2; }
        else if (cps[c] == 2) {
          //outfile_config<<"22 ";
          spinz -= 2; }
        else if (cps[c] == 3) {
          //outfile_config<<"T0 ";
        }
        else if (cps[c] == 4) {
          //outfile_config<<"S0 ";
        }
        else {
          //outfile_config<<"N0 ";
        }
      }
    }
  else {
    if (cps[1] == 1) {
      //outfile_config<<"1 ";
      spinz += 1; }
      else if (cps[1] == 2) {
      //outfile_config<<"2 ";
      spinz -= 1; }
    for (int c = 2; c<= N-1; c+=2){
      if (cps[c] == 1){
        //outfile_config<<"11 ";
        spinz += 2;  }
      else if (cps[c] == 2) {
        //outfile_config<<"22 ";
        spinz -= 2; }
      else if (cps[c] == 3) {
        //outfile_config<<"T0 ";
      }
      else if (cps[c] == 4) {
        //outfile_config<<"S0 ";
      }
      else {
        //outfile_config<<"N0 ";
      }
    }
    if (cps[N] == 1) {
      //outfile_config<<"1 ";
      spinz += 1; }
      else if (cps[N] == 2) {
      //outfile_config<<"2 ";
      spinz -= 1; }
  }
      //outfile_config<<" "<<spinz;

    if (tstep%2 == 1){
            for(int n = 1; n <= N; n = n+2)
             {
               auto s1 = sites(n);
               auto s2 = sites(n+1);
               ITensor upup(s1,s2),dndn(s1,s2),updn1(s1,s2),updn2(s1,s2);
               ITensor pupup(s1,prime(s1),s2,prime(s2)),pdndn(s1,prime(s1),s2,prime(s2));
               ITensor pupdn1(s1,prime(s1),s2,prime(s2)),pupdn2(s1,prime(s1),s2,prime(s2));

               upup.set(s1(1),s2(1),1); dndn.set(s1(2),s2(2),1);
               updn1.set(s1(1),s2(2),sqrt(2)/2); updn1.set(s1(2),s2(1),sqrt(2)/2);
               updn2.set(s1(1),s2(2),sqrt(2)/2); updn2.set(s1(2),s2(1),-sqrt(2)/2);

               pupup.set(s1(1),prime(s1)(1),s2(1),prime(s2)(1),1);
               pdndn.set(s1(2),prime(s1)(2),s2(2),prime(s2)(2),1);

               pupdn1.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
               pupdn1.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
               pupdn1.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),0.5);
               pupdn1.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),0.5);

               pupdn2.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
               pupdn2.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
               pupdn2.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),-0.5);
               pupdn2.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),-0.5);
               //PrintData(upup);

               ITensor U(s1),D,V;
                   if (cps[n] == 1){
                     //cout<<"stepn cps1 "<<(n-1)<<endl;
                     svd(upup,U,D,V,{"Cutoff",cutoff});
                   } else if (cps[n]==2){
                     //cout<<"stepn cps2 "<<(n-1)<<endl;
                     svd(dndn,U,D,V,{"Cutoff",cutoff});
                   } else if (cps[n]==3) {
                     //cout<<"stepn cps3 "<<(n-1)<<endl;
                     svd(updn1,U,D,V,{"Cutoff",cutoff});
                   } else if (cps[n]==4) {
                     //cout<<"stepn cps4 "<<(n-1)<<endl;
                     svd(updn2,U,D,V,{"Cutoff",cutoff});
                   } else {
                     cout<<"init wrong"<<endl;
                   }
                 psi.Anc(n) = U;
                 psi.Anc(n+1) = D*V;
               }
        }
        else
        {

          auto s = sites(1);
          ITensor sitesN1(s);
          if(cps[1] == 1) {
            sitesN1.set(s(1),1);
            sitesN1.set(s(2),0);
            psi.setA(1,sitesN1);
          } else if(cps[1] == 2) {
            sitesN1.set(s(1),0);
            sitesN1.set(s(2),1);
            psi.setA(1,sitesN1);
          } else {
            cout<<"head wrong"<<endl;
          }

          s = sites(N);
          ITensor sitesN2(s);
          if(cps[N] == 1) {
            sitesN2.set(s(1),1);
            sitesN2.set(s(2),0);
            psi.setA(N,sitesN2);
          } else if(cps[N] == 2){
            sitesN2.set(s(1),0);
            sitesN2.set(s(2),1);
            psi.setA(N,sitesN2);
          } else {
            cout<<"tail wrong"<<endl;
          }

          for(int n = 2; n <= N-1; n = n+2)
           {
             auto s1 = sites(n);
             auto s2 = sites(n+1);
             ITensor upup(s1,s2),dndn(s1,s2),updn1(s1,s2),updn2(s1,s2);
             ITensor pupup(s1,prime(s1),s2,prime(s2)),pdndn(s1,prime(s1),s2,prime(s2));
             ITensor pupdn1(s1,prime(s1),s2,prime(s2)),pupdn2(s1,prime(s1),s2,prime(s2));

             upup.set(s1(1),s2(1),1); dndn.set(s1(2),s2(2),1);
             updn1.set(s1(1),s2(2),sqrt(2)/2); updn1.set(s1(2),s2(1),sqrt(2)/2);
             updn2.set(s1(1),s2(2),sqrt(2)/2); updn2.set(s1(2),s2(1),-sqrt(2)/2);

             pupup.set(s1(1),prime(s1)(1),s2(1),prime(s2)(1),1);
             pdndn.set(s1(2),prime(s1)(2),s2(2),prime(s2)(2),1);

             pupdn1.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
             pupdn1.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
             pupdn1.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),0.5);
             pupdn1.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),0.5);

             pupdn2.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
             pupdn2.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
             pupdn2.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),-0.5);
             pupdn2.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),-0.5);
             //PrintData(upup);

             ITensor U(s1),D,V;
                 if (cps[n] == 1){
                   //cout<<"stepn cps1 "<<(n-1)<<endl;
                   svd(upup,U,D,V,{"Cutoff",cutoff});
                 } else if (cps[n]==2){
                   //cout<<"stepn cps2 "<<(n-1)<<endl;
                   svd(dndn,U,D,V,{"Cutoff",cutoff});
                 } else if (cps[n]==3) {
                   //cout<<"stepn cps3 "<<(n-1)<<endl;
                   svd(updn1,U,D,V,{"Cutoff",cutoff});
                 } else if (cps[n]==4) {
                   //cout<<"stepn cps4 "<<(n-1)<<endl;
                   svd(updn2,U,D,V,{"Cutoff",cutoff});
                 } else {
                   cout<<"init wrong"<<endl;
                 }
               psi.Anc(n) = U;
               psi.Anc(n+1) = D*V;
             }
        }
     //PrintData(psi);

     for (int node = 1; node <= N; node++){
       cps[node] = 0;
      }

    //=============================initial state===============================
    exactApplyMPO(psi,expH,psi);
    for(int n = 2; n <= nt; ++n)
        {
        fitApplyMPO(psi,expH,psi,args);
        //printfln("step:%2.2d Energy = %.20f",n,psiHphi(psi,H,psi));
        }
    Energy = psiHphi(psi,H,psi);
    printfln("Energy = %.20f",Energy);

    //===============================collapse===================================
    if (tstep%2 == 0) {
        for(int n = 1; n <= N; n = n+2)
         {
          if(n==1)  psi.position(n);
          ITensor wf = psi.Anc(n)*psi.Anc(n+1);
          auto s1 = sites(n);
          auto s2 = sites(n+1);
          ITensor upup(s1,s2),dndn(s1,s2),updn1(s1,s2),updn2(s1,s2);
          ITensor pupup(s1,prime(s1),s2,prime(s2)),pdndn(s1,prime(s1),s2,prime(s2));
          ITensor pupdn1(s1,prime(s1),s2,prime(s2)),pupdn2(s1,prime(s1),s2,prime(s2));

          upup.set(s1(1),s2(1),1); dndn.set(s1(2),s2(2),1);
          updn1.set(s1(1),s2(2),sqrt(2)/2); updn1.set(s1(2),s2(1),sqrt(2)/2);
          updn2.set(s1(1),s2(2),sqrt(2)/2); updn2.set(s1(2),s2(1),-sqrt(2)/2);

          pupup.set(s1(1),prime(s1)(1),s2(1),prime(s2)(1),1);
          pdndn.set(s1(2),prime(s1)(2),s2(2),prime(s2)(2),1);

          pupdn1.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
          pupdn1.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
          pupdn1.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),0.5);
          pupdn1.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),0.5);

          pupdn2.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
          pupdn2.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
          pupdn2.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),-0.5);
          pupdn2.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),-0.5);

          //get_rand = (rand()%(100-0+1)/100.0);
          get_rand = (rand()%(10000-0+1)/10000.01);
          //cout<<"rand: "<<get_rand<<endl;
          //PrintData(wf);
          Psuu = (dag(prime(wf,Site))*pupup*wf).real();
          Psdd = (dag(prime(wf,Site))*pdndn*wf).real();
          Psud1 = (dag(prime(wf,Site))*pupdn1*wf).real();
          Psud2 = (dag(prime(wf,Site))*pupdn2*wf).real();

          //printfln("upup: %.10f dndn = %.10f updn1: %.10f updn2 = %.10f sum = %.10f", Psuu, Psdd, Psud1, Psud2, Psuu+Psdd+Psud1+Psud2);

          //PrintData(psi.Anc(n+2));
          if (get_rand<=Psuu){
            cps[n]=1;
            if (n<N-1) psi.Anc(n+2) = (upup*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psuu);
            //cout<<"here1: "<<Psuu<<endl;
          } else if (Psuu<get_rand && get_rand<=(Psuu+Psdd)) {
            cps[n]=2;
            if (n<N-1) psi.Anc(n+2) = (dndn*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psdd);
            //cout<<"here2: "<<(Psuu+Psdd)<<endl;
          } else if ((Psuu+Psdd)<get_rand && get_rand<=(1 - Psud2)) {
            cps[n]=3;
            if (n<N-1) psi.Anc(n+2) = (updn1*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psud1);
            //cout<<"here3: "<<(1 - Psud2)<<endl;
          } else if ((1 - Psud2)<get_rand && get_rand<=1) {
            cps[n]=4;
            if (n<N-1) psi.Anc(n+2) = (updn2*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psud2);
            //cout<<"here4: "<<1<<endl;
          } else {
              cout<<"something wrong"<<endl;
            }
          //cout<<"cps: "<<cps[n]<<endl;
        }
      }  else {

            psi.position(1);
            ITensor wf = psi.Anc(1);
            auto s = sites(1);
            ITensor psim1(s),psim2(s);
            psim1.set(s(1),1);psim1.set(s(2),0);
            psim2.set(s(1),0);psim2.set(s(2),1);
            //get_rand = (rand()%(100-0+1)/100.0);
            get_rand = (rand()%(10000-0+1)/10000.01);
            auto Pszu = (dag(prime(wf,Site))*sites.op("projUp",1)*wf).real();
            auto Pszd = (dag(prime(wf,Site))*sites.op("projDn",1)*wf).real();
              if (Pszu >= get_rand)
              { cps[1]=1;
                psi.Anc(1+1) = (psim1*psi.Anc(1)*psi.Anc(1+1))/sqrt(Pszu); }
              else
              { cps[1]=2;
                psi.Anc(1+1) = (psim2*psi.Anc(1)*psi.Anc(1+1))/sqrt(Pszd);
              }
              //cout<<"cps: "<<cps[1]<<endl;
              //printfln(" zup: %.10f zdn = %.10f sum = %.10f", Pszu, Pszd, Pszu+Pszd);

        for(int n = 2; n <= N-1; n = n+2)
         {
          //cout<<"shift2: "<<n<<endl;
          //if(n==1)  psi.position(n);
          ITensor wf = psi.Anc(n)*psi.Anc(n+1);
          auto s1 = sites(n);
          auto s2 = sites(n+1);
          ITensor upup(s1,s2),dndn(s1,s2),updn1(s1,s2),updn2(s1,s2);
          ITensor pupup(s1,prime(s1),s2,prime(s2)),pdndn(s1,prime(s1),s2,prime(s2));
          ITensor pupdn1(s1,prime(s1),s2,prime(s2)),pupdn2(s1,prime(s1),s2,prime(s2));

          upup.set(s1(1),s2(1),1); dndn.set(s1(2),s2(2),1);
          updn1.set(s1(1),s2(2),sqrt(2)/2); updn1.set(s1(2),s2(1),sqrt(2)/2);
          updn2.set(s1(1),s2(2),sqrt(2)/2); updn2.set(s1(2),s2(1),-sqrt(2)/2);

          pupup.set(s1(1),prime(s1)(1),s2(1),prime(s2)(1),1);
          pdndn.set(s1(2),prime(s1)(2),s2(2),prime(s2)(2),1);

          pupdn1.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
          pupdn1.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
          pupdn1.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),0.5);
          pupdn1.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),0.5);

          pupdn2.set(s1(1),s2(2),prime(s1)(1),prime(s2)(2),0.5);
          pupdn2.set(s1(2),s2(1),prime(s1)(2),prime(s2)(1),0.5);
          pupdn2.set(s1(1),s2(2),prime(s1)(2),prime(s2)(1),-0.5);
          pupdn2.set(s1(2),s2(1),prime(s1)(1),prime(s2)(2),-0.5);

          //get_rand = (rand()%(100-0+1)/100.0);
          get_rand = (rand()%(10000-0+1)/10000.01);
          //cout<<"rand: "<<get_rand<<endl;
          //PrintData(wf);
          Psuu = (dag(prime(wf,Site))*pupup*wf).real();
          Psdd = (dag(prime(wf,Site))*pdndn*wf).real();
          Psud1 = (dag(prime(wf,Site))*pupdn1*wf).real();
          Psud2 = (dag(prime(wf,Site))*pupdn2*wf).real();

          //printfln("upup: %.10f dndn = %.10f updn1: %.10f updn2 = %.10f sum = %.10f", Psuu, Psdd, Psud1, Psud2, Psuu+Psdd+Psud1+Psud2);

          //PrintData(psi.Anc(n+2));
          if (get_rand<=Psuu){
            cps[n]=1;
            if (n<N-1) psi.Anc(n+2) = (upup*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psuu);
          } else if (Psuu<get_rand && get_rand<=(Psuu+Psdd)) {
            cps[n]=2;
            if (n<N-1) psi.Anc(n+2) = (dndn*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psdd);
          } else if ((Psuu+Psdd)<get_rand && get_rand<=(1 - Psud2)) {
            cps[n]=3;
            if (n<N-1) psi.Anc(n+2) = (updn1*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psud1);
          } else if ((1 - Psud2)<get_rand && get_rand<=1) {
            cps[n]=4;
            if (n<N-1) psi.Anc(n+2) = (updn2*psi.Anc(n)*psi.Anc(n+1)*psi.Anc(n+2))/sqrt(Psud2);
          } else {
            cout<<"something wrong"<<endl;
            }
          //cout<<"cps: "<<cps[n]<<endl;
        }

        //cout<<"shift2_tail: "<<N<<endl;
        wf = psi.Anc(N);
        //cout<<"real rand: "<<rand()%(100-0+1)<<endl;
        get_rand = (rand()%(10000-0+1)/10000.01);
        Pszu = (dag(prime(wf,Site))*sites.op("projUp",N)*wf).real();
        Pszd = (dag(prime(wf,Site))*sites.op("projDn",N)*wf).real();
        //cout<<"tail_rand: "<<get_rand<<endl;
        if (Pszu >= get_rand)
          { cps[N]=1;
          }
          else
          { cps[N]=2;
          }
          //printfln(" zup: %.10f zdn = %.10f sum = %.10f", Pszu, Pszd, Pszu+Pszd);
          //cout<<"cps: "<<cps[N]<<endl;

      }

      //outfile_config<<" "<<Energy<<" "<<tstep<<endl;
    if (tstep>heatstep)
    {
      Energysum += Energy;
      ii += 1;
      binning_sum[test-1] += Energy;
      //cout<< "Energyave "<<ii<<" "<<Energysum/ii<<endl;
    }

    outfile_thermal_2<<tstep<<" "<<Energy<<" "<<Energysum/ii<<endl;

  }

}

  for (int n =0; n<= test_num-1; n++) {
    binning_ave[n] = binning_sum[n]/((thermal_step - heatstep));
    cout<<n<<" binning_ave: "<<binning_ave[n]<<endl;
  }

  Finite_energy = Energysum/(test_num*(thermal_step - heatstep));
  cout<<"Energysum "<<ii<<" "<<Finite_energy<<endl;
  //cout<<" temperature: "<<20*0.1<< " E: "<<Finite_energy <<endl;
  cout<<" temperature: "<<nt*0.1<< " E: "<<Finite_energy <<endl;

  auto binning = 0.0;
  for (int n1 = 0;n1<=test_num-1;n1++) {
    binning += pow((binning_ave[n1]-Finite_energy),2);
    cout<<"n1: "<<n1<<" "<<binning<<" "<<pow((binning_ave[n1]-Finite_energy),2)<<endl;
  }
  variance = sqrt(binning/(test_num-1));
  cout<<"binning_variace: "<<binning<<" "<<variance<<endl;
  outfile_FiniteE<<nt*0.1<<" "<<Finite_energy<<" "<<variance<<endl;

}


  outfile_thermal_2.close();
  outfile_thermal_8.close();
  outfile_ensemble.close();
  outfile_FiniteE.close();
  outfile_config.close();

return 0;
}
