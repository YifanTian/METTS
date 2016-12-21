#include "itensor/all.h"
#include <stdio.h>      /* printf */
#include <math.h>       /* tanh, log */

using namespace itensor;
using namespace std;


Real c1(int m, Real M){
  auto x = m/M;
  auto c1y = 0.5*(1-tanh((x-0.5)/(x*(1-x))));
  return c1y;
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
    auto N = sites.N();
    auto M = 10.0;
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
        cout<<"i: "<<b<<" "<<ti<<endl;
          H.at(b) = ti*sites.op("Sz",b)*sites.op("Sz",b+1)
                     + 0.5*sites.op("S+",b)*sites.op("S-",b+1)
                     + 0.5*sites.op("S-",b)*sites.op("S+",b+1);

        }
    return H;
    }

int
main(int argc, char* argv[])
    {
    ofstream outfile_bond,sbs_coeff;
    outfile_bond.open("energybond_sbs.txt"),sbs_coeff.open("coefficient.txt");

    int N = 60;

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = SpinHalf(N); //make a chain of N spin 1/2's
    //auto sites = SpinOne(N); //make a chain of N spin 1's

    //
    // Use the AutoMPO feature to create the
    // next-neighbor Heisenberg model
    //
    auto ampo = AutoMPO(sites);
    int M = 10;
    auto ti = 0.0;
    for(int j = 1; j < N; ++j)
        {

          if (j>= 1 && j<=M) {
            ti = c1(M-j,M); }
          else if (j>M && j<=N-M){
            ti = 1;
          } else if(j>N-M && j<N){
            ti = c1(j-N+M,M);
          }
        //cout<<"i: "<<b<<" "<<ti<<endl;
        ampo += ti*0.5,"S+",j,"S-",j+1;
        ampo += ti*0.5,"S-",j,"S+",j+1;
        ampo += ti,    "Sz",j,"Sz",j+1;
        }
    auto H = MPO(ampo);
    auto effH = makeH(sites);
    auto effH0 = makeH0(sites);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i)
        {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
        }
    auto psi = MPS(state);

    //
    // overlap calculates matrix elements of MPO's with respect to MPS's
    // overlap(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", overlap(psi,H,psi) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep.
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(10);
    sweeps.maxm() = 10,20,50,50,75,75,100,100,200,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f", overlap(psi,H,psi) );

// SvN
    psi.position(N/2);
  	auto wf = psi.A(N/2)*psi.A(N/2+1);
  	auto U = psi.A(N/2);
  	ITensor S,V;
  	auto spectrum = svd(wf,U,S,V);

  	Real SvN = 0.0;

  	for(auto p:spectrum.eigs())
  	{
  		if(p>1E-12) SvN += -p*log(p);
  	}
  	printfln("Across bond b=%d, SvN = %.10f",N/2,SvN);

    // energybond
    auto Energy1 = 0.0;
    auto energybond = 0.0;
    int center = int(N/2);
    for(int i = 1;i < N; i++) {
    psi.position(i);
    auto phi = psi.A(i)*psi.A(i+1);
    auto bra = prime(phi,Site);
    //energybond = (bra*effH.at(i)*phi).real();
    energybond = (bra*effH0.at(i)*phi).real();
    Energy1 += energybond; //measurement of energy <S.S>
    outfile_bond<<energybond<<" ";
    }
    outfile_bond<<endl;


    for(int i = 1;i<=60;i++){
      auto cm = c1(i,40.0);
      sbs_coeff<<cm<<" ";
    }
    sbs_coeff<<endl;


    outfile_bond.close();
    sbs_coeff.close();

    return 0;
    }
