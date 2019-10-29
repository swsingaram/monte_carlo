BeginPackage["MonteCarlo`"]

Needs["GraphUtilities`"];


MCEnsembleOneD::usage = "MCEnsembleOneD[m,n,no,trials,step,e] runs the Metropolis Algorithm on m 1D lattices each with n sites.  The total # of proteins in the system is no.  trials is the number of trials and step is the number of trials between two samples, e is the lateral energy"

MCBindAndMoveGraph::usage = "MCBindAndMoveGraph[adjmatfile,occupiedsites,linkedlist,numsteps,numbreaths,epsilon] randomly selects a vertex to bind a protein, then moves the protein around numsteps number of times following the Metropolis Algorithm.  epsilon is the interaction energy, occupiedsites is the list of occupied sites, and adjmatfile is the RAG ourput. the linked list keeps track of next-nearest neighbors and numbreaths is the number of times the program samples this interaction. It contains Degree Matrix, so make sure you take this out first. Remove the brackets and commas too."

MCBindAndMoveAndTransferT3::usage = "Monte Carlo simulation on protein transfer of two nucleo-protein complexes.  The large RNA is represented by a T=3 graph of dimers and the
small RNA is represented by a T=2 graph of dimers.  The list of edges decribing these graphs must be loaded in the run script because the function takes graphs as inputs.  Basically everything is same as MCBindAndMove with the addition of a protein transfer step."

MCBindAndMoveAndTransferT2::usage = "Monte Carlo simulation on protein transfer of two nucleo-protein complexes.  The large RNA dimer is represented by a T=2 graph of dimers and the
small RNA monomer is represented by a T=2 graph of dimers.  The list of edges decribing these graphs must be loaded in the run script because the function takes graphs as inputs.  Basically everything is same as MCBindAndMove with the addition of a protein transfer step."

MC2DJump::usage = "MC governed jumps between two lattices.  The inputs for this function usually comes from MC2DTranslate.  Input-RNAmat1, RNAmat2, Energy, epsilon."

MC2DAddParticle::usage = "Add a particle to the 2D lattice and return new energy based on nearest-neighbor interactions.  Input-RNAmat, Energy, epsilon."

MC2DTranslate::usage = "MC simulation on 2D lattice. Translations only.   Prepare RNA matrix with defective sites.  0-unoccupied, 1-occupied, 2-defective.  As input, give the matrix rnamat, a list of the occupied sites (i.e. {{2,2},{1,4},etc}), the number of steps in the MC, the nearest-neighbor interaction energy epsilon, and the total energy"

MCMoveGraph::usage = "Monte Carlo simulation for protein transfer of two nucleo-protein complexes.  The large RNA is represented by a T=3/T=2 graph of dimers. This function performs the translations onthe capsid/deformed capsid.  INPUT-graph of capsid/deformed capsid, list of occupied vertices, number of translation moves, the nearest-neighbor energy, the total energy."

MCJumpGraph::usage = "Monte Carlo simulation for protein transfer of two nucleo-protein complexes.  The large RNA is typically represented by a T=3 graph of dimers and the small RNA is represented by a T=2 graph of dimers. This function performs a monte-carlo governed jump between the two capsids/deformed capsids.  INPUT-graph of capsid1/deformed capsid1, list of occupied vertices1, graph of capsid1/deformed capsid2,list of vertices2, the nearest-neighbor energy, the total energy."

MCAddGraph::usage = "Monte Carlo simulation for adding protein to graph.  RNA is typically represnted by a graph of dimers.  This function randomly picks an unoccupied site and fills, then calcuates the change in energy.  INPUT - graph of capsids, list of occupied sites, nearest-neighbor interaction energy, total energy."

MCHopGraph::usage = "Monte Carlo simulation for letting protein hop to other unoccupied sites on the RNA.  RNA is represented as a graph of dimers.  This function randomly picks an occupied site and an unoccupied site, then switches them according to the MC algorithm.  INPUT -graph of capsids, list of occupied sites, nearest-neighbor interaction energy, total energy."

RNAGraph::usage = "RNAGraph means we remove vertices from the graph so that the number of lattice sites is the same as the number of binding sites.  INPUT--graph, list deleted vertices, list of removed edges, number of vertices to delete.  OUTPUT--graph, list of deleted vertices, list of removed edges"

LowerpH::usage = "lowering the pH means we add vertices to the graph.  INPUT--graph, list of deleted vertices, list of removed edges, number of vertices to add.  OUTPUT--graph, list of deleted vertices, list of removed edges"

StaggerGraph::usage = "staggers the vertices of a graph.  For example, if the vertices are {1,2,3,4 } this function would return {5,6,7,8}, so that one can combine the two graphs using the GraphUnion mathematica function.  INPUT--graph.  OUTPUT-staggered graph."

RNAAssociate::usage = "only used to model RNA association via CP dimers.  INPUT--graph1 (union of two graphs, where edges are staggard), occsites, valency1,valency2,epsilon, energy.  OUTPUT-graph union"

RNAIslandAssociate::usage = "brings islands from multicomponent graphs together.  Just looks at the occupied sites.  INPUT--graph,occsutes, associated edges, valency,epsilon, energy.  OUTPUT-graph" 




Begin["`Private`"]

MCEnsembleOneD[m_,n_,no_,trials_,step_,e_] := (OccupiedSites = {}; While[Less[Length[OccupiedSites], no],OccupiedSites = Union[OccupiedSites,DeleteDuplicates[RandomSample[Range[Times[n,m]], no]]];UnOccupiedSites = Complement[Range[Times[n,m]], OccupiedSites];Config = ConstantArray[0,Times[n,m]];Energy = 0;Do[Config[[OccupiedSites[[i]]]] = 1;If[And[Unequal[Mod[OccupiedSites[[i]],n],1], Equal[Config[[OccupiedSites[[i]] - 1]], 1]],Energy = Energy - e;];,{i,1,no}];Do[TempConfig = Config;TempOcc = OccupiedSites;  TempUnOcc = UnOccupiedSites; randI1 = RandomInteger[{1,no}];  randOcc = OccupiedSites[[randI1]];  randI2 = RandomInteger[{1, Times[m,n] - no}];  randUnOcc = UnOccupiedSites[[randI2]];  TempConfig[[randOcc]] = 0;  TempConfig[[randUnOcc]] = 1; TempOcc[[randI1]] = randUnOcc;  TempUnOcc[[randI2]] = randOcc;  Dele = 0;  nni = 0;  nnf = 0;    If[And[Unequal[Mod[randOcc,n],0],Equal[Config[[randOcc + 1]], 1]],  nni = nni + 1;,  nni = nni;  ];    If[And[Unequal[Mod[randOcc,n],1],Equal[Config[[randOcc - 1]], 1]],      nni = nni + 1;,      nni = nni;      ];    If[And[Unequal[Mod[randUnOcc,n],0],Equal[TempConfig[[randUnOcc + 1]], 1]],      nnf = nnf + 1;,      nnf = nnf;     ];  If[And[Unequal[Mod[randUnOcc,n],1],Equal[TempConfig[[randUnOcc - 1]], 1]],  nnf = nnf + 1;,  nnf = nnf;  ];Delnn = nnf - nni;Dele = Times[Delnn,Times[-1,e]];  If[LessEqual[Dele,0],    Energy = Energy + Dele;   Config = TempConfig;    OccupiedSites = TempOcc;    UnOccupiedSites = TempUnOcc;,    r = RandomReal[];    If[Less[r,Exp[Times[-1,Dele]]],       Energy = Energy + Dele;       Config = TempConfig;       OccupiedSites = TempOcc;       UnOccupiedSites = TempUnOcc;,       Energy = Energy;       Config = Config;       OccupiedSites = OccupiedSites;       UnOccupiedSites = UnOccupiedSites;      ];     ];If[Equal[Mod[numtrials,step],0], Do[    WriteConfig = Part[Config,((i - 1)*n + 1);;(i*n)]; NumProt = Count[WriteConfig,1]; PutAppend[WriteConfig,StringJoin[ToString[Divide[no,m]],"configurations.dat"]]; PutAppend[NumProt,StringJoin[ToString[Divide[no,m]],"numprotein.dat"]]; ,{i,1,m}];  PutAppend[Energy,StringJoin[ToString[Divide[no,m]],"energy.dat"]]; ];,{numtrials,1,trials}];];)

RNAAssociate[dimergraph_,occsites_,associatededges_,stagger_,valency_,epsilon_,energy_] := (

DimerGraph=dimergraph;

OccSites=occsites;

Stagger=stagger;

Valency=valency;

AssocEdges=associatededges;

Energy=energy;





(*only one vertex per RNA can be associated to another*)
VertAssocE= VertexList[Graph[AssocEdges]];


OccSites1=Intersection[OccSites,Range[1,Stagger]];

OccSites2=Intersection[OccSites,Range[Stagger+1,2*Stagger]];

If[OccSites1 == {} || OccSites2 == {},
       (*cannot associate*)
        Return[{DimerGraph,OccSites,AssocEdges,Valency,Energy}];
  ];

(*the remaining occupied sites on each RNA that can associate*)

OccSites1=Complement[OccSites1,Intersection[OccSites1,VertAssocE]];

OccSites2=Complement[OccSites2,Intersection[OccSites2,VertAssocE]];




If[OccSites1 != {} && OccSites2 != {},
       
   	randOcc1=OccSites1[[RandomInteger[{1,Length[OccSites1]}]]];
        
	randOcc2=OccSites2[[RandomInteger[{1,Length[OccSites2]}]]];
        NNOccSites1=Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc1,1],{randOcc1}],OccSites1];
	NNUnOccSites1=Complement[Complement[NeighborhoodVertices[DimerGraph,randOcc1,1],{randOcc1}],NNOccSites1];
	NNOccSites2=Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc2,1],{randOcc2}],OccSites2];
	NNUnOccSites2=Complement[Complement[NeighborhoodVertices[DimerGraph,randOcc2,1],{randOcc2}],NNOccSites2];

       	If[Length[NNOccSites1] < Valency && Length[NNOccSites2] < Valency,
    	 DimerGraph=EdgeAdd[DimerGraph,UndirectedEdge[randOcc1,randOcc2]];
     	AssocEdges=Append[AssocEdges,UndirectedEdge[randOcc1,randOcc2]];
     	Energy=Energy+epsilon;
		If[VertexDegree[DimerGraph,randOcc1] > Valency,
                    (*delete the unoccupied edge*)
                    randUnOcc1=NNUnOccSites1[[RandomInteger[{1,Length[NNUnOccSites1]}]]];
                    DimerGraph=EdgeDelete[DimerGraph,UndirectedEdge[randOcc1,randUnOcc1]];    
		];
                If[VertexDegree[DimerGraph,randOcc2] > Valency,
                    (*delete the unoccupied edge*)
                    randUnOcc2=NNUnOccSites2[[RandomInteger[{1,Length[NNUnOccSites2]}]]];
                    DimerGraph=EdgeDelete[DimerGraph,UndirectedEdge[randOcc2,randUnOcc2]];    
		];
        ];
 
    ];

(*now randomly choose an edge to break*)

If[AssocEdges=={},
 Return[{DimerGraph,OccSites,AssocEdges,Valency,Energy}];
  ];    

randAssocEdge=AssocEdges[[RandomInteger[{1,Length[AssocEdges]}]]];

(*pick a random real number between 0 & 1*)
r= RandomReal[];
If[r<Exp[epsilon],
            (*accept move*)
             Energy=Energy-epsilon;
             AssocEdges=Complement[AssocEdges,{randAssocEdge}];
             DimerGraph=EdgeDelete[DimerGraph,randAssocEdge];
 ];             

Return[{DimerGraph,OccSites,AssocEdges,Valency,Energy}];

)

RNAIslandAssociate[dimergraph_,occsites_,associatededges_,valency_,epsilon_,energy_] := (

DimerGraph=dimergraph;

OccSites=occsites;

AssocEdges=associatededges;

Valency=valency;

Energy=energy;

(*maybe you should restrict which guys can associate*)

VertAssocE=VertexList[Graph[AssocEdges]];

OccSitesSubGraph=Subgraph[DimerGraph,Complement[OccSites,VertAssocE]];

ConnectedCompOcc=ConnectedComponents[OccSitesSubGraph];

If[Length[ConnectedCompOcc] == 1,
       Return[{DimerGraph,OccSites,AssocEdges,Energy}];
  ];

(*Now pick occupied two occupied sites from different islands*)
r1=RandomInteger[{1,Length[ConnectedCompOcc]}];
r2=RandomInteger[{1,Length[ConnectedCompOcc]}];

While[ r2 == r1,

 r2=RandomInteger[{1,Length[ConnectedCompOcc]}];

     ];

island1=ConnectedCompOcc[[r1]];
island2=ConnectedCompOcc[[r2]];

randOcc1=island1[[RandomInteger[{1,Length[island1]}]]];
randOcc2=island2[[RandomInteger[{1,Length[island2]}]]];
(*Put the RNA associate part*)
NNOccSites1=Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc1,1],{randOcc1}],OccSites];
	NNUnOccSites1=Complement[Complement[NeighborhoodVertices[DimerGraph,randOcc1,1],{randOcc1}],NNOccSites1];
	NNOccSites2=Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc2,1],{randOcc2}],OccSites];
	NNUnOccSites2=Complement[Complement[NeighborhoodVertices[DimerGraph,randOcc2,1],{randOcc2}],NNOccSites2];

       	If[Length[NNOccSites1] < Valency && Length[NNOccSites2] < Valency,
    	 DimerGraph=EdgeAdd[DimerGraph,UndirectedEdge[randOcc1,randOcc2]];
     	AssocEdges=Append[AssocEdges,UndirectedEdge[randOcc1,randOcc2]];
     	Energy=Energy+epsilon;
		If[VertexDegree[DimerGraph,randOcc1] > Valency,
                    (*delete the unoccupied edge*)
                    randUnOcc1=NNUnOccSites1[[RandomInteger[{1,Length[NNUnOccSites1]}]]];
                    DimerGraph=EdgeDelete[DimerGraph,UndirectedEdge[randOcc1,randUnOcc1]];    
		];
                If[VertexDegree[DimerGraph,randOcc2] > Valency,
                    (*delete the unoccupied edge*)
                    randUnOcc2=NNUnOccSites2[[RandomInteger[{1,Length[NNUnOccSites2]}]]];
                    DimerGraph=EdgeDelete[DimerGraph,UndirectedEdge[randOcc2,randUnOcc2]];    
		];
        ];

(*now randomly choose an edge to break*)

If[AssocEdges=={},
 Return[{DimerGraph,OccSites,AssocEdges,Energy}];
  ];    

randAssocEdge=AssocEdges[[RandomInteger[{1,Length[AssocEdges]}]]];

(*pick a random real number between 0 & 1*)
r= RandomReal[];
If[r<Exp[epsilon],
            (*accept move*)
             Energy=Energy-epsilon;
             AssocEdges=Complement[AssocEdges,{randAssocEdge}];
             DimerGraph=EdgeDelete[DimerGraph,randAssocEdge];
 ]; 

 Return[{DimerGraph,OccSites,AssocEdges,Energy}];

)


StaggerGraph[dimergraph_,stagger_]:= (

DimerGraph=dimergraph;

Elist=EdgeList[DimerGraph];

NumVertices=stagger;

NewElist={};


(*for every edge, get vertices, stagger them both *)
Do[
(*Vertices of the edge*)
VerticesElist=VertexList[Graph[{i}]];

(*Now stagger each vertex*)
   Do[
      VerticesElist[[j]]=VerticesElist[[j]]+NumVertices;
      ,{j,1,2}];
NewElist=Append[NewElist,UndirectedEdge[VerticesElist[[1]],VerticesElist[[2]]]];
,{i,Elist}];

DimerGraphStagger=Graph[NewElist];

Return[DimerGraphStagger];


)


RNAGraph[dimergraph_,vdel_,edel_,numvdel_]:= (

DimerGraph=dimergraph;

OccSites=occsites;

Vdel=vdel;

Edel=edel;

NumVdel=numvdel;

(*we will randomly remove vertices, but not the ones which are occupied*)

RemoveableVertices=VertexList[DimerGraph];

(*randomly choose numvdel to delete*)
posns={};
While[Length[posns] < numvdel,
   posns = Union[posns,DeleteDuplicates[RandomSample[Range[Length[RemoveableVertices]],numvdel]]];
      ];

(*for each vertex we delete, we need to keep track of the edges*)
Edel1={};
Vdel1={};
Do[
(*find neighborhood vertices*)
nhvert=Complement[NeighborhoodVertices[DimerGraph,RemoveableVertices[[i]],1],{RemoveableVertices[[i]]}];
nhedges={};
	Do[
	nhedges=Append[nhedges,UndirectedEdge[RemoveableVertices[[i]],j]];
	,{j,nhvert}];

Vdel1=Append[Vdel1,RemoveableVertices[[i]]];
Edel1=Append[Edel1,nhedges];

,{i,posns}];

DimerGraph=VertexDelete[DimerGraph,Vdel1];

Vdel=Join[Vdel,Vdel1];
Edel=Join[Edel,Edel1];

Return[{DimerGraph,Vdel,Edel}];

)

LowerpH[dimergraph_,vdel_,edel_,numvadd_]:= (

DimerGraph=dimergraph;
(*deleted vertices*)
Vdel=vdel;

(*deleted edges*)
Edel=edel;

(*number of vertices to add*)
NumVadd=numvadd;

addposns={};

While[Length[addposns] < NumVadd,
addposns = Union[addposns,DeleteDuplicates[RandomSample[Range[Length[Vdel]],NumVadd]]];
     ]; 
(*We'll make the EdgeList*)

edgelist={};

Do[
edgelist=Append[edgelist,Edel[[i]]];
DimerGraph=VertexAdd[DimerGraph,Vdel[[i]]];

,{i,addposns}];

(*take out the edge repeats*)
edgelist=Flatten[edgelist];

edgelist=Union[Sort /@ edgelist];


DimerGraph=EdgeAdd[DimerGraph,edgelist];
Vdel=Delete[Vdel,Partition[addposns,1]];
Edel=Delete[Edel,Partition[addposns,1]];

Return[{DimerGraph,Vdel,Edel}];

)


MCAddGraph[dimergraph_,occsites_,epsilon_,energy_]:= (

DimerGraph=dimergraph;

OccSites=occsites;

Energy=energy;

UnOccSites=Complement[VertexList[DimerGraph],OccSites];


If[Length[UnOccSites] == 0,
     Return[{DimerGraph,OccSites,Energy}];
  ];

(*pick a random unoccupied site to occupy*)
randUnOcc = UnOccSites[[RandomInteger[{1,Length[UnOccSites]}]]];
(*occupied sites next to it/energy*)
OccNN = Intersection[Complement[NeighborhoodVertices[DimerGraph,randUnOcc,1],{randUnOcc}],OccSites];
DelE= Length[OccNN]*epsilon;

Energy=Energy+DelE;
UnOccSites=Complement[UnOccSites,{randUnOcc}];
OccSites=Append[OccSites,randUnOcc];

Return[{DimerGraph,OccSites,Energy}];

)


MCHopGraph[dimergraph_,occsites_,epsilon_,energy_] := (

DimerGraph=dimergraph;

OccSites=occsites;

Energy=energy;

UnOccSites=Complement[VertexList[DimerGraph],OccSites];

If[Length[OccSites] == 0 || Length[UnOccSites] == 0,
     Return[{DimerGraph,OccSites,Energy}];
  ];

(*randomly pick an occupied site*)
randOcc = OccSites[[RandomInteger[{1,Length[OccSites]}]]];
(* pick a random unoccupied site*)
randUnOcc = UnOccSites[[RandomInteger[{1,Length[UnOccSites]}]]];
(*calculate initial energy*)
E1 = Length[Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc,1],{randOcc}],OccSites]]*epsilon;
(*calculate final energy*)
E2 = Length[Intersection[NeighborhoodVertices[DimerGraph,randUnOcc,1],Complement[OccSites,{randOcc}]]]*epsilon;

DelE = E2-E1;
 (*Now Determine whether we accept or reject the move and update OccSites, randOcc, and UnOccSites*)
   If[DelE <= 0,
          (*Accept Move!*)
          OccSites = Complement[Append[OccSites,randUnOcc],{randOcc}];
	  UnOccSites = Complement[Append[UnOccSites,randOcc],{randUnOcc}];
          Energy = Energy + DelE;
                    ,
       
         (*DelEnergy > 0, so accept boltzman probability*)
         r = RandomReal[];
         If[r<Exp[-(DelE)],
                   (*Accept Move*)
                   OccSites = Complement[Append[OccSites,randUnOcc],{randOcc}];                   
                   UnOccSites = Complement[Append[UnOccSites,randOcc],{randUnOcc}];
	           Energy = Energy + DelE;	
                            ,
                
                   (*Reject Move*)  
                   OccSites = OccSites;
                   UnOccSites = UnOccSites;
                   Energy = Energy;
                 
           ];
                    
    ];

Return[{DimerGraph,OccSites,Energy}];


)



MCMoveGraph[dimergraph_,occsites_,epsilon_,energy_] := (

DimerGraph=dimergraph;

OccSites=occsites;

Energy=energy;

UnOccSites=Complement[VertexList[DimerGraph],OccSites];


If[Length[OccSites] == 0 || Length[UnOccSites] == 0,
     Return[{DimerGraph,OccSites,Energy}];
  ];


(*randomly pick an occupied site*)
randOcc = OccSites[[RandomInteger[{1,Length[OccSites]}]]];
(* unoccupied sites next to it *)
UnOccNN = Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc,1],{randOcc}],UnOccSites];
If[UnOccNN == {},
     (*Can't Move*)
     OccSites = OccSites;
     UnOccSites = UnOccSites;
     Energy = Energy;
     ,

    (*calculate the initial energy.  I think you only need to calculate it for the local environment*)
    InitialEnergy = Length[Complement[Complement[NeighborhoodVertices[DimerGraph,randOcc,1],{randOcc}],UnOccNN]]*epsilon;

   (*pick a random unoccupied nearest-neighbor*)
   randUnOccNN = UnOccNN[[RandomInteger[{1,Length[UnOccNN]}]]];
  
   (*figure out the occupied nearest-neighbors of the randomly chosen nearest-neighbors NOT including the randOcc*)
   OccNNNN = Intersection[Complement[NeighborhoodVertices[DimerGraph,randUnOccNN,1],{randUnOccNN}],Complement[OccSites,{randOcc}]];
   (*If we moved the occupied site to the unoccupied site, calculate the final energy*)
   FinalEnergy = Length[OccNNNN]*epsilon;
   DelEnergy = FinalEnergy - InitialEnergy;
   (*Now Determine whether we accept or reject the move and update OccSites, randOcc, and UnOccSites*)
   If[DelEnergy <= 0,
          (*Accept Move!*)
          OccSites = Complement[Append[OccSites,randUnOccNN],{randOcc}];
	  UnOccSites = Complement[Append[UnOccSites,randOcc],{randUnOccNN}];
          Energy = Energy + DelEnergy;
                    ,
       
         (*DelEnergy > 0, so accept boltzman probability*)
         r = RandomReal[];
         If[r<Exp[-(DelEnergy)],
                   (*Accept Move*)
                   OccSites = Complement[Append[OccSites,randUnOccNN],{randOcc}];                   
                   UnOccSites = Complement[Append[UnOccSites,randOcc],{randUnOccNN}];
	           Energy = Energy + DelEnergy;	
                            ,
                
                   (*Reject Move*)  
                   OccSites = OccSites;
                   UnOccSites = UnOccSites;
                   Energy = Energy;
                 
           ];
                    
    ];
  ];




Return[{DimerGraph,OccSites,Energy}];

)


MCJumpGraph[dimergraph1_,occsites1_,dimergraph2_,occsites2_,epsilon_,energy_]:=(

DimerGraph1=dimergraph1;

OccSites1=occsites1;

Energy=energy;

UnOccSites1=Complement[VertexList[DimerGraph1],OccSites1];

DimerGraph2=dimergraph2;

OccSites2=occsites2;

UnOccSites2=Complement[VertexList[DimerGraph2],OccSites2];

OccSites=Join[OccSites1,OccSites2];

If[OccSites != {},
    (*pick a random occupied site*)
    randOcc=OccSites[[RandomInteger[{1,Length[OccSites]}]]];
   
	If[MemberQ[OccSites1,randOcc],
                (*randOcc is on dimergraph1*)
 	         randOcc1=randOcc;

	(*first we'll choose a random OCCUPIED site on dimergraph1 and an UNOCCUPIED site on dimergraph2*)

	(*If[Length[OccSites1] == 0 || Length[UnOccSites2] == 0,
     	 Return[{DimerGraph1,OccSites1,DimerGraph2,OccSites2,DimerGraph2,Energy}];
	  ];*)
	If[Length[OccSites1] > 0 && Length[UnOccSites2] > 0,
	(*Pick a random occupied site dimergraph1 and a random unoccupied site from dimergraph2*)
     	(*randOcc1=OccSites1[[RandomInteger[{1,Length[OccSites1]}]]];*)
     	randUnOcc2 = UnOccSites2[[RandomInteger[{1,Length[UnOccSites2]}]]];
     
     
     	(*calculate initial energy*)
     	OccNN = Intersection[Complement[NeighborhoodVertices[DimerGraph1,randOcc1,1],{randOcc1}],OccSites1];
   
      
     	InitialEnergy = Length[OccNN]*epsilon;
     

	(*calculate final energy*)
    	OccNN = Intersection[Complement[NeighborhoodVertices[DimerGraph2,randUnOcc2,1],{randUnOcc2}],OccSites2];
    
    
    	FinalEnergy = Length[OccNN]*epsilon;

    	DelEnergy = FinalEnergy - InitialEnergy;
    
   

     		If[DelEnergy <= 0,
          	(*Accept Move*)
          	OccSites2 = Append[OccSites2,randUnOcc2];
         	 UnOccSites2 = Complement[UnOccSites2,{randUnOcc2}];
          	OccSites1 = Complement[OccSites1,{randOcc1}];
          	UnOccSites1 = Append[UnOccSites1,randOcc1];

          	Energy = Energy + DelEnergy;
          	,
          	(*DelEnergy > 0, accept with boltzman probability*)
          	r = RandomReal[];
          		If[r<Exp[-(DelEnergy)],
               		(*Accept Move*)
               		OccSites2 = Append[OccSites2,randUnOcc2];
               		UnOccSites2 = Complement[UnOccSites2,{randUnOcc2}];
              		OccSites1 = Complement[OccSites1,{randOcc1}];
               		UnOccSites1 = Append[UnOccSites1,randOcc1];

              		 Energy = Energy + DelEnergy;
             		 ,
             		 (*Reject Move*)
             		 Energy = Energy;
          		 ];
        	];
	];
 
        ,
         (*randOcc is  on dimergraph2*)
          randOcc2 = randOcc;
       (*Pick a random OCCUPIED site from DimerGraph2 and a random UNOCCUPIED site from DimerGraph1*)
	
	(*If[Length[OccSites2] == 0 || Length[UnOccSites1] == 0,
      	Return[{DimerGraph1,OccSites1,DimerGraph2,OccSites2,DimerGraph2,Energy}];
 	 ];*)

	If[Length[OccSites2] > 0 && Length[UnOccSites1] > 0,
       (*randOcc2 = OccSites2[[RandomInteger[{1,Length[OccSites2]}]]];*)
       	randUnOcc1 = UnOccSites1[[RandomInteger[{1,Length[UnOccSites1]}]]];

       (*calculate initial energy*)
       OccNN = Intersection[Complement[NeighborhoodVertices[DimerGraph2,randOcc2,1],{randOcc2}],OccSites2];
       InitialEnergy = Length[OccNN]*epsilon;
       (*calculate final energy*)
       OccNN = Intersection[Complement[NeighborhoodVertices[DimerGraph1,randUnOcc1,1],{randUnOcc1}],OccSites1];
      
       FinalEnergy = Length[OccNN]*epsilon;

       DelEnergy = FinalEnergy - InitialEnergy;
     
     

      		If[DelEnergy <= 0,
          	(*Accept Move*)
          	OccSites1 = Append[OccSites1,randUnOcc1];
          	UnOccSites1 = Complement[UnOccSites1,{randUnOcc1}];
         	 OccSites2 = Complement[OccSites2,{randOcc2}];
          	UnOccSites2 = Append[UnOccSites2,randOcc2];

          	Energy = Energy + DelEnergy;
         	 ,
         	 (*DelEnergy > 0, accept with boltzman probability*)
         	 r = RandomReal[];
         		 If[r<Exp[-(DelEnergy)],
             	 	(*Accept Move*)
             		 OccSites1 = Append[OccSites1,randUnOcc1];
             	 	UnOccSites1 = Complement[UnOccSites1,{randUnOcc1}];
             		 OccSites2 = Complement[OccSites2,{randOcc2}];
             		 UnOccSites2 = Append[UnOccSites2,randOcc2];

             		 Energy = Energy + DelEnergy;
              		,
              		(*Reject Move*)
             		 Energy = Energy;
           		];
      		 ];

	 ];

     ];

];
Return[{DimerGraph1,OccSites1,DimerGraph2,OccSites2,Energy}];


)


MCBindAndMoveGraph[adjmat_,occupiedsites_, linkedlist_,numsteps_,numbreaths_,epsilon_] := (
(*adjmatfile is the RAG ourput.  It contains Degree Matrix, so make sure you take this out first. Remove the brackets and commas too.*)
(*AdjMatList = Import[adjmatfile,"Table"];*)
(*Make Tree graph from adjacency matrix*)
AdjMatList = adjmat;
AdjGraph = AdjacencyGraph[AdjMatList];
(*check to see if tree graph is connected*)
(*If[TreeGraphQ[AdjGraph] == False, Return["Tree graph representation of secondary structure is not connected."];];*)

(*
no = Length[occupiedsites];

(*File names*)
occsitesfilename = StringJoin[ToString[no], "_OccSitesconfig.dat"];
*)

(*vertices of the unoccupied sites*)
UnoccupiedSites = Complement[VertexList[AdjGraph],occupiedsites];

(*Pick a random unoccupied site, this site is now occupied*)
RandomOcc = UnoccupiedSites[[RandomInteger[{1,Length[UnoccupiedSites]}]]];
OccupiedSites = Append[occupiedsites,RandomOcc];
UnoccupiedSites = Complement[UnoccupiedSites,{RandomOcc}];
TaggedOcc = RandomOcc;
(*Metropolois Algorithm*)
Energy=0;
LinkedVertices=linkedlist;
Do[

(*figure out the unoccupied nearest-neighbors*)
UnOccNN = Intersection[Complement[NeighborhoodVertices[AdjGraph,TaggedOcc,1],{TaggedOcc}],UnoccupiedSites];
If[UnOccNN == {},
                (*Can't Move*)
        OccupiedSites = OccupiedSites;
        UnoccupiedSites = UnoccupiedSites;
        TaggedOcc = TaggedOcc;
        Energy = Energy;
,

	(*calculate the initial energy.  I think you only need to calculate it for the local environment*)
	InitialEnergy = Length[Complement[Complement[NeighborhoodVertices[AdjGraph,TaggedOcc,1],{TaggedOcc}],UnOccNN]]*epsilon;

	(*pick a random unoccupied nearest-neighbor*)
	TaggedUnOcc = UnOccNN[[RandomInteger[{1,Length[UnOccNN]}]]];


	(*figure out the occupied nearest-neighbors of the randomly chosen nearest neighbor NOT including the RandOcc *)
	OccNNNN = Intersection[Complement[NeighborhoodVertices[AdjGraph,TaggedUnOcc,1],{TaggedUnOcc}],Complement[OccupiedSites,{TaggedOcc}]];
	(*If we moved the occupied site to the unoccupied site, calculate the final energy*)
	FinalEnergy = Length[OccNNNN]*epsilon;

	DelEnergy = FinalEnergy - InitialEnergy;

(*Now determine whether we accept or reject the move and update OccupiedSites, TagOcc, and UnOccupiedSites*)

	If[DelEnergy <= 0,
		(*Accept Move!*)
		OccupiedSites = Complement[Append[OccupiedSites,TaggedUnOcc],{TaggedOcc}];
		UnoccupiedSites = Complement[Append[UnoccupiedSites,TaggedOcc],{TaggedUnOcc}];
		TaggedOcc = TaggedUnOcc;
		Energy = Energy + DelEnergy;
	,
		(*DelEnergy > 0, so accept with boltzman probability *)
		r = RandomReal[];
		If[r<Exp[-(DelEnergy)],
   			(*Accept Move*)
			OccupiedSites = Complement[Append[OccupiedSites,TaggedUnOcc],{TaggedOcc}];
			UnoccupiedSites = Complement[Append[UnoccupiedSites,TaggedOcc],{TaggedUnOcc}];
			TaggedOcc = TaggedUnOcc;
			Energy = Energy + DelEnergy;
		,
			(*Reject Move*)
			OccupiedSites = OccupiedSites;
			UnoccupiedSites = UnoccupiedSites;
			TaggedOcc = TaggedOcc;
			Energy = Energy;
		];
	];
	(*
	********************
	********************
	ALLOW PROTEINS
	TWO LINKS AWAY TO
	INTERACT.  THIS MODELS
	RNA BREATHING.
	********************
	********************
	*)
	If[Mod[i,numbreaths] == 0,
   		(*Allow a breath every once in a while*)
                 (*
                 *******************
                 *******************
                 ALLOW TAGGEDOCC
                 TO ESCAPE
                 *******************
                 *******************
                 *)
                 (*Using LinkedVertices, a select a random pair*)
                 
                If[LinkedVertices == {},
                      (*Can't Move*)
                       AdjGraph = AdjGraph;,
                      
		       ri=RandomInteger[{1,Length[LinkedVertices]}];
                       RandPair = LinkedVertices[[ri]];
                      (*Now allow an escape if the random number is less than the Bolzmann probability*)
                   r = RandomReal[];
				If[ r <= Exp[(epsilon)],
                                           (*Accept Move*)
                                           (*Break Bond*)
                                           AdjMatList[[RandPair[[1]],RandPair[[2]]]]=0;
                                           AdjMatList[[RandPair[[2]],RandPair[[1]]]]=0;
					   AdjGraph=AdjacencyGraph[AdjMatList];
					   LinkedVertices=Delete[LinkedVertices,ri];
                                  ,
         					
					(*Reject Move*)
                                        AdjGraph = AdjGraph;
                                  ];

                   ];

  		(*Pick a random occupied site, figure out the occupied sites TWO links away*)
		RandOcc1=OccupiedSites[[RandomInteger[{1,Length[OccupiedSites]}]]];
  		OccSitesTwoLinks = Intersection[Complement[NeighborhoodVertices[AdjGraph,RandOcc1,2],NeighborhoodVertices[AdjGraph,RandOcc1,1]],OccupiedSites];
   
		If[OccSitesTwoLinks == {},
          		(*Can't Move*)
           		AdjGraph = AdjGraph;,

			RandOcc2 = OccSitesTwoLinks[[RandomInteger[{1,Length[OccSitesTwoLinks]}]]];
			(*Determine if we allow the two proteins to interact by treating
			the three-vertex subgraph as an ideal polymer and accepting the move
			based on end-to-end ideal polymer statistics.*)
			(*So I'll accept the move if the end-to-end distance is less than or equal to 2 nm.  I've assumed a duplex is 1.5 nm (5bp and 0.3 nm/bp).*)
   			(*The probability ends up being 0.46*)
			r = RandomReal[];
    				If[r <= 0.46,
        				(*Accept Move *)
        				(*Add link between two vertices*)
         				AdjMatList[[RandOcc1,RandOcc2]]=1;
         				AdjMatList[[RandOcc2,RandOcc1]]=1;
         				AdjGraph=AdjacencyGraph[AdjMatList];
					LinkedVertices=AppendTo[LinkedVertices,{RandOcc1,RandOcc2}];
        			,

        				(*Reject Move*)
        				AdjGraph = AdjGraph;



     				 ];

 		 ];
	];



];

,{i,1,numsteps}];

 
Join[{OccupiedSites},{AdjMatList},{LinkedVertices}]

)

MCBindAndMoveAndTransferT3[t3dimergraph_,t2dimergraph_,numprott3_,numprott2_,transfertry_,nummoves_,reportstep_,epsilon_]:=(


(* put this part in your run script as well as T2
(*Import each of the graphs*)
(*T=3 graph of dimers*)
t3dimerlist=Import[t3dimergraphfile];
t3dimer=Graph[t3dimerlist,Rule[DirectedEdges,False]];

*)

t3dimer=t3dimergraph;
t2dimer=t2dimergraph;


(*uncomment if you want 3D graphics
t3dimerVcoor=GraphCoordinates3D[t3dimerlist];*)

(* put this is in run script
(*T=2 graph of dimers*)
t2dimerlist=Import[t2dimergraphfile];
t2dimer=Graph[t2dimerlist,Rule[DirectedEdges,False]];
*)

(*uncomment if you want 3D graphics
t2dimerVcoor=GraphCoordinates3D[t2dimerlist];*)

(*First we'll allow proteins on each molecule to move around a bit before the molecules
collide *)

(*occupied sites of T=3*)
T3OccSites={};
(*randomly populate sites on T=3*)
While[Length[T3OccSites] < numprott3,
T3OccSites = Union[T3OccSites,DeleteDuplicates[RandomSample[Range[90],numprott3]]]];

(*occupied sites of T=2*)
T2OccSites={};
(*randomly populate sites on T=2*)
While[Length[T2OccSites] < numprott2,
T2OccSites = Union[T2OccSites,DeleteDuplicates[RandomSample[Range[60],numprott2]]]];

(*unoccupied sites of T=3*)
T3UnOccSites=Complement[Range[90],T3OccSites];

(*unoccupied sites of T=2*)
T2UnOccSites=Complement[Range[60],T2OccSites];

Energy=0;

Do[
(************************************
*************************************
MOVES ON T=3
*************************************
************************************)
If[Length[T3OccSites] == 0 || Length[T2OccSites] == 0,
     Return["Finished!"]];
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
(*(* uncomment this section if you want 3D graphics *)
***************************************************
****3D GRAPHICS**********************************
(*First, we need to make a list of coordinates of the occupied sites*)
T3command={};
Do[
coordinate=t3dimerVcoor[[i]];
T3command=Append[T3command,Red];
T3command=Append[T3command,Sphere[coordinate,0.1]];
,{i,T3OccSites}];

(*Now do the same for T=2*)
T2command={};
Do[
coordinate=t2dimerVcoor[[i]];
T2command=Append[T2command,Red];
T2command=Append[T2command,Sphere[coordinate,0.1]];
,{i,T2OccSites}];
Print["T=3 command"];
Print[T3command];
Print["T=2 command"];
Print[T2command];

(*Export[StringJoin[ToString[movenum],"_t3graph.pdf"],GraphPlot3D[t3dimerlist,Rule[VertexRenderingFunction,(T3command &)],Rule[PlotStyle,{Black}]]];*)
(*Export[StringJoin[ToString[movenum],"_t2graph.pdf"],GraphPlot3D[t2dimerlist,Rule[VertexRenderingFunction,(T2command &)],Rule[PlotStyle,{Black}]]];*)
************************************************************
************************************************************
*)
  
 Export[StringJoin[ToString[movenum],"_t3graph.pdf"],HighlightGraph[t3dimer,Subgraph[t3dimer,T3OccSites]]];
Export[StringJoin[ToString[movenum],"_t2graph.pdf"],HighlightGraph[t2dimer,Subgraph[t2dimer,T2OccSites]]]; 
];
,{movenum,1,nummoves}];

)

MCBindAndMoveAndTransferT2[t3dimergraph_,t2dimergraph_,numprott3_,numprott2_,transfertry_,nummoves_,reportstep_,epsilon_]:=(

(* put this part in run script
(*Import each of the graphs*)
(*T=3 graph of dimers*)

t3dimer=Graph[Import[t3dimergraphfile]];
*)
(*T=2 graph of dimers*)
(*put this in run script
t2dimer=Graph[Import[t2dimergraphfile]];
*)
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
Export[StringJoin[ToString[movenum],"_t3graph.pdf"],HighlightGraph[t3dimer,Subgraph[t3dimer,T3OccSites]]];
Export[StringJoin[ToString[movenum],"_t2graph.pdf"],HighlightGraph[t2dimer,Subgraph[t2dimer,T2OccSites]]];
  ];


,{movenum,1,nummoves}];

)

MC2DJump[rnamat1_,rnamat2_,epsilonrna1_,epsilonrna2_,energy_] := (

(*make local variables*)
RNA1mat=rnamat1;
RNA2mat=rnamat2;
Energy=energy;

NumRowsRNA1=Dimensions[RNA1mat][[1]];
NumColumnsRNA1=Dimensions[RNA1mat][[2]];

NumRowsRNA2=Dimensions[RNA2mat][[1]];
NumColumnsRNA2=Dimensions[RNA2mat][[2]];

OccupiedSitesRNA1=Position[RNA1mat,1];
UnOccupiedSitesRNA2=Position[RNA2mat,0];

If[ OccupiedSitesRNA1 == {} || UnOccupiedSitesRNA2 == {},
       (*can't move*)
       Return[{RNA1mat,RNA2mat,Energy}];
  ];

(*pick random occupied site on RNA1 and unoccupiedsite RNA2*)
RandomOccRNA1=OccupiedSitesRNA1[[RandomInteger[{1,Length[OccupiedSitesRNA1]}]]];
RandomUnOccRNA2=UnOccupiedSitesRNA2[[RandomInteger[{1,Length[UnOccupiedSitesRNA2]}]]];





(*Calculate initial energy*)
(*we must determine the occupied and unoccupied nearest-neighbors*)
NNOccSitesRNA1=0;


(*RNA*)

	Do[
		If[RNA1mat[[Mod[RandomOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomOccRNA1[[2]]+j,NumColumnsRNA1,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomOccRNA1[[2]]+j,NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA1mat[[Mod[RandomOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomOccRNA1[[2]],NumColumnsRNA1,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomOccRNA1[[2]],NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

InitialEnergy=NNOccSitesRNA1*epsilonrna1;




NNOccSitesRNA2=0;
	Do[
		If[RNA2mat[[Mod[RandomUnOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]]+j,NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomUnOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]]+j,NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA2mat[[Mod[RandomUnOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]],NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomUnOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]],NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

FinalEnergy=NNOccSitesRNA2*epsilonrna2;

DelE=FinalEnergy-InitialEnergy;

(*MC time!*)
If[DelE <= 0,
         (*accept move*)
	 RNA2mat[[RandomUnOccRNA2[[1]],RandomUnOccRNA2[[2]]]]=1;
         RNA1mat[[RandomOccRNA1[[1]],RandomOccRNA1[[2]]]]=0;
	Energy=Energy+DelE;
	,
	(*accept move with boltzman probability*)
       	r = RandomReal[];
	If[r<Exp[-(DelE)],
             (*accept move*)	
	 	RNA2mat[[RandomUnOccRNA2[[1]],RandomUnOccRNA2[[2]]]]=1;
        	RNA1mat[[RandomOccRNA1[[1]],RandomOccRNA1[[2]]]]=0;
		Energy=Energy+DelE;
		,
		(*reject move*)
		Energy=Energy;
	];

];


OccupiedSitesRNA2=Position[RNA2mat,1];
UnOccupiedSitesRNA1=Position[RNA1mat,0];

(*pick random occupied site on RNA2 and unoccupied site on RNA1*)
If[ OccupiedSitesRNA2 == {} || UnOccupiedSitesRNA1 == {},
       (*can't move*)
       Return[{RNA1mat,RNA2mat,Energy}];
  ];

(*pick random occupied site on RNA2 and unoccupiedsite RNA1*)
RandomOccRNA2=OccupiedSitesRNA2[[RandomInteger[{1,Length[OccupiedSitesRNA2]}]]];
RandomUnOccRNA1=UnOccupiedSitesRNA1[[RandomInteger[{1,Length[UnOccupiedSitesRNA1]}]]];



(*Calculate initial energy*)
(*we must determine the occupied and unoccupied nearest-neighbors*)
NNOccSitesRNA2=0;


(*RNA*)

	Do[
		If[RNA2mat[[Mod[RandomOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomOccRNA2[[2]]+j,NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomOccRNA2[[2]]+j,NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA2mat[[Mod[RandomOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomOccRNA2[[2]],NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomOccRNA2[[2]],NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

InitialEnergy=NNOccSitesRNA2*epsilonrna2;



NNOccSitesRNA1=0;
	Do[
		If[RNA1mat[[Mod[RandomUnOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]]+j,NumColumnsRNA1,1]]] != 2,
		(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomUnOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]]+j,NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA1mat[[Mod[RandomUnOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]],NumColumnsRNA1,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomUnOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]],NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

FinalEnergy=NNOccSitesRNA1*epsilonrna1;

DelE=FinalEnergy-InitialEnergy;

(*MC time!*)
If[DelE <= 0,
         (*accept move*)
	 RNA2mat[[RandomOccRNA2[[1]],RandomOccRNA2[[2]]]]=0;
         RNA1mat[[RandomUnOccRNA1[[1]],RandomUnOccRNA1[[2]]]]=1;
	Energy=Energy+DelE;
	,
	(*accept move with boltzman probability*)
       	r = RandomReal[];
	If[r<Exp[-(DelE)],
             (*accept move*)	
	 	RNA2mat[[RandomOccRNA2[[1]],RandomOccRNA2[[2]]]]=0;
        	RNA1mat[[RandomUnOccRNA1[[1]],RandomUnOccRNA1[[2]]]]=1;
		Energy=Energy+DelE;
		,
		(*reject move*)
		Energy=Energy;
	];

];
 

Join[{RNA1mat},{RNA2mat},{Energy}]

)








MC2DAddParticle[rnamat_,epsilon_,energy_] := (
(*make local copies of variables*)
RNAmat=rnamat;
Energy=energy;
NumRows=Dimensions[RNAmat][[1]];
NumColumns=Dimensions[RNAmat][[2]];

(*Find unoccupiedsites*)

UnOccupiedSites=Position[RNAmat,0];

If[UnOccupiedSites == {},
   Return[{RNAmat,Energy}];
  ];

RandomUnOcc=UnOccupiedSites[[RandomInteger[{1,Length[UnOccupiedSites]}]]];

DelE=0;
(*now calculate the number of occupied nearest-neighbors*)
NNOccSites=0;

	Do[
		If[RNAmat[[Mod[RandomUnOcc[[1]],NumRows,1],Mod[RandomUnOcc[[2]]+j,NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomUnOcc[[1]],NumRows,1],Mod[RandomUnOcc[[2]]+j,NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = NNOccSites+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNAmat[[Mod[RandomUnOcc[[1]]+i,NumRows,1],Mod[RandomUnOcc[[2]],NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomUnOcc[[1]]+i,NumRows,1],Mod[RandomUnOcc[[2]],NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = NNOccSites+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

DelE=NNOccSites*epsilon;
Energy=Energy+DelE;
(*update matrix*)
RNAmat[[RandomUnOcc[[1]],RandomUnOcc[[2]]]]=1;

Join[{RNAmat},{Energy}]

)












MC2DTranslate[rnamat_, numsteps_,epsilon_,energy_] := (


(*make local variables of RNAmat,occupiedsitesRNA,occupiedsitesRNA2*)
RNAmat=rnamat;
NumRows=Dimensions[RNAmat][[1]];
NumColumns=Dimensions[RNAmat][[2]];
Energy=energy;

OccupiedSites=Position[RNAmat,1];

UnOccupiedSites=Position[RNAmat,0];

If[OccupiedSites == {} || UnOccupiedSites == {},
   Return[{RNAmat,Energy}];];

(*LOOP TIME BABY*)

Do[
(*pick random occupied sites on RNA *)

RandomOcc=OccupiedSites[[RandomInteger[{1,Length[OccupiedSites]}]]];

(*we must determine the occupied and unoccupied nearest-neighbors*)
NNOccSites={};
NNUnOccSites={};

(*RNA*)

	Do[
		If[RNAmat[[Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = Append[NNOccSites,{Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		NNUnOccSites = Append[NNUnOccSites,{Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]}];	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNAmat[[Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = Append[NNOccSites,{Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		NNUnOccSites = Append[NNUnOccSites,{Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]}];	
			];

		];

	,{i,-1,1,2}];
(*initial energy*)
InitialEnergy=Length[NNOccSites]*epsilon;
If[NNUnOccSites == {},
       (*can't move*)
	FinalEnergy=InitialEnergy;
	OccNNNN={};
			,
	
	(*not empty, so make a virtual move*)	
         (*pick random unoccupied nearest-neighbor*)
	  RandomNNUnOccSites=NNUnOccSites[[RandomInteger[{1,Length[NNUnOccSites]}]]];
         (*figure out the occupied nearest-neighbors of the randomly chosen unoccupied nearest neighbor NOT including RandOcc*)
	OccNNNN={};
	
	Do[
		If[RNAmat[[Mod[RandomNNUnOccSites[[1]],NumRows,1],Mod[RandomNNUnOccSites[[2]]+j,NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomNNUnOccSites[[1]],NumRows,1],Mod[RandomNNUnOccSites[[2]]+j,NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					OccNNNN = Append[OccNNNN,{Mod[RandomNNUnOccSites[[1]],NumRows,1],Mod[RandomNNUnOccSites[[2]]+j,NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		(*we don't need to keep track of them*)
					Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNAmat[[Mod[RandomNNUnOccSites[[1]]+i,NumRows,1],Mod[RandomNNUnOccSites[[2]],NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomNNUnOccSites[[1]]+i,NumRows,1],Mod[RandomNNUnOccSites[[2]],NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					OccNNNN = Append[OccNNNN,{Mod[RandomNNUnOccSites[[1]]+i,NumRows,1],Mod[RandomNNUnOccSites[[2]],NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		(*we don't need to keep track of them*)	
					Energy=Energy;
			];

		];

	,{i,-1,1,2}];



	(*but this procedure will have included RandomOcc*)

        
	FinalEnergy=Length[Complement[OccNNNN,{RandomOcc}]]*epsilon;

  ];

(*Now calculate change in Energy*)
DelE=FinalEnergy-InitialEnergy;

(*Now make monte-carlo move*)
If[DelE <= 0,
   (*accept move*)
  	  If[NNUnOccSites != {},
	    (*we can move*)
	    Energy=Energy+DelE;
	    (*swap occupied site for the randomly chosen unoccupied site*)
	     OccupiedSites=Complement[Append[OccupiedSites,RandomNNUnOccSites],{RandomOcc}];
	     UnOccupiedSites=Complement[Append[UnOccupiedSites,RandomOcc],{RandomNNUnOccSites}];
	(*update matrix*)
	    RNAmat[[RandomNNUnOccSites[[1]],RandomNNUnOccSites[[2]]]]=1;
	    RNAmat[[RandomOcc[[1]],RandomOcc[[2]]]]=0;
		,
	     (*can't move*)
	     Energy=Energy;
	     OccupiedSites=OccupiedSites;
	     UnOccupiedSites=UnOccupiedSites;
	    ];
	,

	(*DelE is > 0 so accept with boltzman proabability*)
	r = RandomReal[];
	If[r<Exp[-(DelE)],
		(*accept move*)
		Energy=Energy+DelE;
		OccupiedSites=Complement[Append[OccupiedSites,RandomNNUnOccSites],{RandomOcc}];
	        UnOccupiedSites=Complement[Append[UnOccupiedSites,RandomOcc],{RandomNNUnOccSites}];
	    (*update matrix*)
	    RNAmat[[RandomNNUnOccSites[[1]],RandomNNUnOccSites[[2]]]]=1;
	    RNAmat[[RandomOcc[[1]],RandomOcc[[2]]]]=0;
		,
		(*reject move*)
		Energy=Energy;
		OccupiedSites=OccupiedSites;
		UnOccupiedSites=UnOccupiedSites;


          ];

 ];

,{i,1,numsteps}];

Join[{RNAmat},{Energy}]

)





End[ ]

EndPackage[ ]
