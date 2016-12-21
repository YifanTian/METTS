#include "itensor/all.h"

using namespace std;
using namespace itensor;
using std::vector;


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

int main()
{
int N = 50,upn,dnn,spinz;
//int N = 10,upn,dnn,spinz;
auto sites = SpinHalf(N);
auto state = InitState(sites);
Real cutoff = 1E-12;
Real get_rand, Psuu,Psdd,Psud1,Psud2,Pszu,Pszd,Psxu,Psxd,Energysum = 0.0,Finite_energy = 0.0,variance = 0.0,Energy,rsz,rszsum;
//int thermal_step = 400,heatstep = 50,test_num = 5,ii;
int thermal_step = 1049,heatstep = 50,ii,snum = 10;
Real cps[N+1];
Real Energybond[N],qnbond[N+1],qnbondsum[N+1],Energybondsum[N];
//Real binning_ave[test_num],binning_sum[test_num];
Real binning_ave[snum],binning_sum[snum];
Real Etstep[thermal_step];
for(int i = 1;i<=thermal_step;i++) Etstep[i] = 0.0;

ofstream outfile_ensemble,outfile_FiniteE,outfile_config,outfile_C,outfile_Sn,outfile_SS;
ofstream outfile_bond,outfile_qnbond;
outfile_FiniteE.open("Metts_FE.txt");
outfile_config.open("configuration.txt");
outfile_bond.open("bondmix_sbs.txt");
outfile_qnbond.open("qnbond_mix_sbs.txt");

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

auto effH = makeH(sites);
auto effH0 = makeH0(sites);
auto tau = 0.05;
auto expH = toExpH<ITensor>(ampo,tau);

auto args = Args("Cutoff=",1E-9,"Maxm=",3000);
auto ttotal = 1.0;
//auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
auto nt = 10;

for (int ent = 20; ent <= 80; ent=ent+10) {
  nt = ent;

  ii = 0;
  Energysum = 0.0;
  for(int i = 0;i<=snum-1;i++) {
    binning_ave[i] = 0.0;
    binning_sum[i] = 0.0;
  }

  for(int i = 1;i < N; i++){
    Energybondsum[i] = 0.0; //measurement of energy <S.S>
  }

  for(int i = 1;i <= N; i++){
    qnbondsum[i] = 0.0; //measurement of energy <S.S>
  }

  rszsum = 0.0;


//for (int test = 1; test<=test_num; test += 1) {

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
      if (cps[c] == 1) {
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
      //cout<<"step: "<<tstep<<" "<<spinz<<endl;

    if (tstep%3 == 1){
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
        else if (tstep%3 == 2)
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
        } else {
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
    //printfln("Energy = %.20f",Energy);

    auto Energy1 = 0.0;
    int center = int(N/2);
    for(int i = center;i < center + 1; i++) {
    psi.position(i);
    auto phi = psi.A(i)*psi.A(i+1);
    auto bra = prime(phi,Site);
    Energy1 += (bra*effH0.at(i)*phi).real(); //measurement of energy <S.S>
    }

    //Energy = 0.0;
    for(int i = 1;i < N; i++){
      psi.position(i);
      auto phi = psi.A(i)*psi.A(i+1);
      auto bra = prime(phi,Site);
      Energybond[i] = (bra*effH0.at(i)*phi).real(); //measurement of energy <S.S>
      //Energy+= Energybond[i];
    }

    rsz = 0.0;
    for(int i = 1;i<=N;i++) {
      psi.position(i);
      auto phi = psi.A(i);
      auto bra = prime(phi,Site);
      qnbond[i] = pow((bra*sites.op("Sz",i)*phi).real(),2); //measurement of energy <S.S>
      rsz += (bra*sites.op("Sz",i)*phi).real();
    }

    //cout<<"rsz: "<<rsz<<endl;
      for(int i = 1;i < N; i++){
          //cout<<Energybond[i]<<" ";
          }


    //===============================collapse===================================
    if (tstep%3 == 0) {
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
      }  else if (tstep%3 == 1) {

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
        }

          wf = psi.Anc(N);
          //cout<<"real rand: "<<rand()%(100-0+1)<<endl;
          get_rand = (rand()%(10000-0+1)/10000.01);
          Pszu = (dag(prime(wf,Site))*sites.op("projUp",N)*wf).real();
          Pszd = (dag(prime(wf,Site))*sites.op("projDn",N)*wf).real();
          //cout<<"tail_rand: "<<get_rand<<endl;
          if (Pszu >= get_rand)
            { cps[N]=1; }
            else
            { cps[N]=2; }

        } else {

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
            //get_rand = (rand()%(100-0+1)/100.0);
            get_rand = (rand()%(10000-0+1)/10000.01);
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

      //outfile_config<<" "<<Energy<<" "<<tstep<<endl;
      if (tstep>=heatstep)
      {
        //Energysum += Energy;
        Energysum += Energy1;
        ii += 1;
        int binnum = int(((tstep+50)/100)-1);
        binning_sum[binnum]+= Energy1;
        for(int i = 1;i < N; i++){
          Energybondsum[i] += Energybond[i];
        }
        for(int i = 1;i <= N; i++){
          qnbondsum[i] += qnbond[i];
        }
        rszsum += rsz;
      }
    //outfile_thermal_2<<tstep<<" "<<Energy<<" "<<Energysum/ii<<endl;

  }

  for(int i = 1;i < N; i++){
    Energybondsum[i] = Energybondsum[i]/(thermal_step-heatstep+1);
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
  cout<<" temperature: "<<nt*0.1<< " E: "<<Finite_energy <<endl;

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
    }

  outfile_FiniteE.close();
  outfile_config.close();
  outfile_bond.close();
  outfile_qnbond.close();

return 0;
}
