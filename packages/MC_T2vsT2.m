BeginPackage["test`"]

Needs["GraphUtilities`"];

MCBindAndMoveAndTransfer::usage = "Monte Carlo simulation on protein transfer of two nucleo-protein complexes.  The large RNA is represented by a T=3 graph of dimers and the
small RNA is represented by a T=2 graph of dimers.  The list of edges decribing these graphs must be loaded.  Basically everything is same as MCBindAndMove with the addition of a protein transfer step."

Begin["`Private`"]

MCBindAndMoveAndTransfer[t3dimergraphfile_,t2dimergraphfile_,numprott3_,numprott2_,transfertry_,nummoves_,reportstep_,epsilon_]:=(

(*Import each of the graphs*)
(*T=3 graph of dimers*)

t3dimer=Graph[Import[t3dimergraphfile]];

(*T=2 graph of dimers*)

t2dimer=Graph[Import[t2dimergraphfile]];

(*First we'll allow proteins on each molecule to move around a bit before the molecules
collide *)

(*occupied sites of T=3*)
T3OccSites={};
(*randomly populate sites on T=3*)
While[Length[T3OccSites] < numprott3,
T3OccSites = Union[T3OccSites,DeleteDuplicates[RandomSample[Range[60],numprott3]]]];

(*occupied sites of T=2*)
T2OccSites={};
(*randomly populate sites on T=2*)
While[Length[T2OccSites] < numprott2,
T2OccSites = Union[T2OccSites,DeleteDuplicates[RandomSample[Range[60],numprott2]]]];

(*unoccupied sites of T=3*)
T3UnOccSites=Complement[Range[60],T3OccSites];

(*unoccupied sites of T=2*)
T2UnOccSites=Complement[Range[60],T2OccSites];

Energy=0;

Do[
If[Length[T3OccSites] == 0 || Length[T2OccSites] == 0,
     Return["Finished"];
  ];
(************************************
*************************************
MOVES ON T=3
*************************************
************************************)
(*randomly pick an occupied site on T=3*)
randOccT3 = T3OccSites[[RandomInteger[{1,Length[T3OccSites]}]]];
(* unoccupied sites next to it *)
T3UnOccNN = Intersection[Complement[NeighborhoodVertices[t3dimer,randOccT3,1],{randOccT3}],T3UnOccSites];
If[T3UnOccNN == {},
     (*Can't Move*)
     T3OccSites = T3OccSites;
     T3UnOccSites = T3UnOccSites;
     Energy = Energy;
     ,

    (*calculate the initial energy.  I think you only need to calculate it for the local environment*)
    InitialEnergy = Length[Complement[Complement[NeighborhoodVertices[t3dimer,randOccT3,1],{randOccT3}],T3UnOccNN]]*epsilon;

   (*pick a random unoccupied nearest-neighbor*)
   T3randUnOccNN = T3UnOccNN[[RandomInteger[{1,Length[T3UnOccNN]}]]];
  
   (*figure out the occupied nearest-neighbors of the randomly chosen nearest-neighbors NOT including the randOccT3*)
   T3OccNNNN = Intersection[Complement[NeighborhoodVertices[t3dimer,T3randUnOccNN,1],{T3randUnOccNN}],Complement[T3OccSites,{randOccT3}]];
   (*If we moved the occupied site to the unoccupied site, calculate the final energy*)
   FinalEnergy = Length[T3OccNNNN]*epsilon;
   DelEnergy = FinalEnergy - InitialEnergy;
   (*Now Determine whether we accept or reject the move and update T3OccSites, randOccT3, and T3UnOccSites*)
   If[DelEnergy <= 0,
          (*Accept Move!*)
          T3OccSites = Complement[Append[T3OccSites,T3randUnOccNN],{randOccT3}];
	  T3UnOccSites = Complement[Append[T3UnOccSites,randOccT3],{T3randUnOccNN}];
          Energy = Energy + DelEnergy;
                    ,
       
         (*DelEnergy > 0, so accept boltzman probability*)
         r = RandomReal[];
         If[r<Exp[-(DelEnergy)],
                   (*Accept Move*)
                   T3OccSites = Complement[Append[T3OccSites,T3randUnOccNN],{randOccT3}];                   
                   T3UnOccSites = Complement[Append[T3UnOccSites,randOccT3],{T3randUnOccNN}];
	           Energy = Energy + DelEnergy;	
                            ,
                
                   (*Reject Move*)  
                   T3OccSites = T3OccSites;
                   T3UnOccSites = T3UnOccSites;
                   Energy = Energy;
                 
           ];
                    
    ];
  ];
(*************************************
**************************************
MOVE ON T=2
**************************************
*************************************)
(*randomly pick an occupied site on T=2*)
randOccT2 = T2OccSites[[RandomInteger[{1,Length[T2OccSites]}]]];
(*unoccupied sites next to the occupied site*)
T2UnOccNN = Intersection[Complement[NeighborhoodVertices[t2dimer,randOccT2,1],{randOccT2}],T2UnOccSites];
If[T2UnOccNN == {},
   (*Can't Move*)
   T2OccSites = T2OccSites;
   T2UnOccSites = T2UnOccSites;
   Energy = Energy
   ,
   
   (*calculate the initial energy.  I think you only need to calculate it for the local environment*)
   InitialEnergy = Length[Complement[Complement[NeighborhoodVertices[t2dimer,randOccT2,1],{randOccT2}],T2UnOccNN]]*epsilon;

   (*pick a random unoccupied nearest-neighbor*)
   T2randUnOccNN = T2UnOccNN[[RandomInteger[{1,Length[T2UnOccNN]}]]];
  
   (*figure out the occupied nearest-neighbors of the randomly chosen nearest-neighbors NOT including the randOccT2*)
   T2OccNNNN = Intersection[Complement[NeighborhoodVertices[t2dimer,T2randUnOccNN,1],{T2randUnOccNN}],Complement[T2OccSites,{randOccT2}]];
   (*If we moved the occupied site to the unoccupied site, calculate the final energy*)
   FinalEnergy = Length[T2OccNNNN]*epsilon;

   DelEnergy = FinalEnergy - InitialEnergy;
   (*Now Determine whether we accept or reject the move and update T2OccSites, randOccT2, and T2UnOccSites*)
        If[DelEnergy <= 0,
                (*Accept Move!*)
                T2OccSites = Complement[Append[T2OccSites,T2randUnOccNN],{randOccT2}];
	
	        T2UnOccSites = Complement[Append[T2UnOccSites,randOccT2],{T2randUnOccNN}];
            
                Energy = Energy + DelEnergy;
                ,
       
               (*DelEnergy > 0, so accept boltzman probability*)
               r = RandomReal[];
               If[r<Exp[-(DelEnergy)],
                       (*Accept Move*)
                       T2OccSites = Complement[Append[T2OccSites,T2randUnOccNN],{randOccT2}];                   
                       T2UnOccSites = Complement[Append[T2UnOccSites,randOccT2],{T2randUnOccNN}];
	               Energy = Energy + DelEnergy;	
                       ,
                
                      (*Reject Move*)  
                      T2OccSites = T2OccSites;
                      T2UnOccSites = T2UnOccSites;
                      Energy = Energy;
                 
                  ];
                    
         ];

  ];

(*******************************
********************************
COLLISION
********************************
*******************************)
If[Mod[movenum,transfertry] == 0,
     (*code for MC protein-transfer*)
     (*Pick a random unoccupied site from T=2 and a random occupied site from T=3*)
     randOccT3=T3OccSites[[RandomInteger[{1,Length[T3OccSites]}]]];
     randUnOccT2 = T2UnOccSites[[RandomInteger[{1,Length[T2UnOccSites]}]]];
     (*calculate initial energy*)
     T3OccNN = Intersection[Complement[NeighborhoodVertices[t3dimer,randOccT3,1],{randOccT3}],T3OccSites];
     
     InitialEnergy = Length[T3OccNN]*epsilon;
     

(*calculate final energy*)
    T2OccNN = Intersection[Complement[NeighborhoodVertices[t2dimer,randUnOccT2,1],{randUnOccT2}],T2OccSites];

    FinalEnergy = Length[T2OccNN]*epsilon;

    DelEnergy = FinalEnergy - InitialEnergy;
      If[DelEnergy <= 0,
          (*Accept Move*)
          T2OccSites = Append[T2OccSites,randUnOccT2];
          T2UnOccSites = Complement[T2UnOccSites,{randUnOccT2}];
          T3OccSites = Complement[T3OccSites,{randOccT3}];
          T3UnOccSites = Append[T3UnOccSites,randOccT3];

          Energy = Energy + DelEnergy;
          ,
          (*DelEnergy > 0, accept with boltzman probability*)
          r = RandomReal[];
          If[r<Exp[-(DelEnergy)],
               (*Accept Move*)
               T2OccSites = Append[T2OccSites,randUnOccT2];
               T2UnOccSites = Complement[T2UnOccSites,{randUnOccT2}];
               T3OccSites = Complement[T3OccSites,{randOccT3}];
               T3UnOccSites = Append[T3UnOccSites,randOccT3];

               Energy = Energy + DelEnergy;
              ,
              (*Reject Move*)
              Energy = Energy;
           ];
        ];

       (*Pick a random occupied sites from T=2 and a random unoccupied sites from T=3*)
       randOccT2 = T2OccSites[[RandomInteger[{1,Length[T2OccSites]}]]];
       randUnOccT3 = T3UnOccSites[[RandomInteger[{1,Length[T3UnOccSites]}]]];

       (*calculate initial energy*)
       T2OccNN = Intersection[Complement[NeighborhoodVertices[t2dimer,randOccT2,1],{randOccT2}],T2OccSites];
       InitialEnergy = Length[T2OccNN]*epsilon;
       (*calculate final energy*)
       T3OccNN = Intersection[Complement[NeighborhoodVertices[t3dimer,randUnOccT3,1],{randUnOccT3}],T3OccSites];
      
       FinalEnergy = Length[T3OccNN]*epsilon;

       DelEnergy = FinalEnergy - InitialEnergy;
      If[DelEnergy <= 0,
          (*Accept Move*)
          T3OccSites = Append[T3OccSites,randUnOccT3];
          T3UnOccSites = Complement[T3UnOccSites,{randUnOccT3}];
          T2OccSites = Complement[T2OccSites,{randOccT2}];
          T2UnOccSites = Append[T2UnOccSites,randOccT2];

          Energy = Energy + DelEnergy;
          ,
          (*DelEnergy > 0, accept with boltzman probability*)
          r = RandomReal[];
          If[r<Exp[-(DelEnergy)],
              (*Accept Move*)
              T3OccSites = Append[T3OccSites,randUnOccT3];
              T3UnOccSites = Complement[T3UnOccSites,{randUnOccT3}];
              T2OccSites = Complement[T2OccSites,{randOccT2}];
              T2UnOccSites = Append[T2UnOccSites,randOccT2];

              Energy = Energy + DelEnergy;
              ,
              (*Reject Move*)
              Energy = Energy;
           ];
       ];
  ];

If[Mod[movenum,reportstep] ==0,
(*Export[StringJoin[ToString[movenum],"_t3graph.pdf"],HighlightGraph[t3dimer,Subgraph[t3dimer,T3OccSites]]];
Export[StringJoin[ToString[movenum],"_t2graph.pdf"],HighlightGraph[t2dimer,Subgraph[t2dimer,T2OccSites]]];*)
Print[movenum];
  ];


,{movenum,1,nummoves}];

prottransfer=numprott3-Length[T3OccSites];
PutAppend[prottransfer,StringJoin[ToString[numprott3],"_transfer.dat"]];

)

End[ ]

EndPackage[ ]
