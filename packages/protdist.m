BeginPackage["ProtDist`"]




zExact::usage = "zExact[Np,M,e] is the exact partition function for a lattice of M sites Np proteins and an interaction energy of e.  A negative value of e is attractive."

ZExact::usage = "ZExact[M,mu,e] is the grand canonical partition function for a lattice of M sites, chemical potential mu, and an interaction energy of e.  A negative value of e is attractive."

PExact::usage = "PExact[n,M,mu,e] computes the probability of observing a lattice of M sites covered with n proteins.  The average coverage is fixed by mu and there is an interaction energy of e."

a::usage = "a[e,f] compute the number of normalized clusters for a coverage fraction f and a interaction energy.  POSITIVE values of e are attractive."

A::usage = "A[f,e,l] computes the number of normalized clusters of length l for a coverage f and an interaction energy e.  POSITIVE values of e are attractive."

 Begin["`Private`"]
 
 


zExact[Np_,M_,e_] := (Return[N[Sum[Binomial[Np-1,Ng-1]*Binomial[M-Np+1,Ng]*Exp[-(Np-Ng)*e],{Ng,1,Np}]]];)

ZExact[M_,mu_,e_] := (Return[Sum[N[Exp[mu*n]]*zExact[n,M,e],{n,1,M}]];)

PExact[n_,M_,mu_,e_] := (Return[N[(Exp[mu*n]*zExact[n,M,e])/ZExact[M,mu,e]]];)

a[e_,f_] := (Return[(Sqrt[1+4*f*(1-f)*(Exp[e]-1)]-1)/(2*(Exp[e]-1))];)

A[f_,e_,l_] := (Return[((a[e,f])^2/f)*(1-(a[e,f]/f))^(l-1)];)

End[]

EndPackage[]
