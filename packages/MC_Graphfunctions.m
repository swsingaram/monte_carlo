BeginPackage["test`"]

Needs["GraphUtilities`"];

MCMoveGraph::usage = "Monte Carlo simulation for protein transfer of two nucleo-protein complexes.  The large RNA is represented by a T=3/T=2 graph of dimers. This function performs the translations onthe capsid/deformed capsid.  INPUT-graph of capsid/deformed capsid, list of occupied vertices, number of translation moves, the nearest-neighbor energy, the total energy."

MCJumpGraph::usage = "Monte Carlo simulation for protein transfer of two nucleo-protein complexes.  The large RNA is typically represented by a T=3 graph of dimers and the small RNA is represented by a T=2 graph of dimers. This function performs a monte-carlo governed jump between the two capsids/deformed capsids.  INPUT-graph of capsid1/deformed capsid1, list of occupied vertices1, graph of capsid1/deformed capsid2,list of vertices2, the nearest-neighbor energy, the total energy."

MCAddGraph::usage = "Monte Carlo simulation for adding protein to graph.  RNA is typically represnted by a graph of dimers.  This function randomly picks an unoccupied site and fills, then calcuates the change in energy.  INPUT - graph of capsids, list of occupied sites, nearest-neighbor interaction energy, total energy."

RaisepH::usage = "raising the pH means we remove vertices from the graph.  INPUT--graph, occsites, list deleted vertices, list of removed edges, number of vertices to delete.  OUTPUT--graph, list of deleted vertices, list of removed edges"

LowerpH::usage = "lowering the pH means we add vertices to the graph.  INPUT--graph, list of deleted vertices, list of removed edges, number of vertices to add.  OUTPUT--graph, list of deleted vertices, list of removed edges"

(*EdgeUppH::usage = "raising the pH means we remove vertices from graph, keeping vertices the same.  INPUT--graph,occsites,list of removed edges, energy, epsilon, number of edges to delete.  OUTPUT-graph, list of removed edges, energy, epsilon"*)

StaggerGraph::usage = "staggers the vertices of a graph.  For example, if the vertices are {1,2,3,4 } this function would return {5,6,7,8}, so that one can combine the two graphs using the GraphUnion mathematica function.  INPUT--graph.  OUTPUT-staggered graph."

RNAAssociate:usage = "only used to model RNA association via CP dimers.  INPUT--graph1 (union of two graphs, where edges are staggard), occsites, valency1,valency2,epsilon, energy.  OUTPUT-graph union"


Begin["`Private`"]

RNAAssociate[dimergraph_,occsites_,associatededges_,tnum_,valency_,epsilon_,energy_] := (

DimerGraph=dimergraph;

OccSites=occsites;

Tnum=tnum;

Valency=valency;

AssocEdges=associatededges;

Energy=energy;


(*number of vertices on a graph is 30 times the T number*)
NumVertices=30*Tnum;

(*only one vertex per RNA can be associated to another*)
VertAssocE= VertexList[Graph[AssocEdges]];


OccSites1=Intersection[OccSites,Range[1,NumVertices]];

OccSites2=Intersection[OccSites,Range[NumVertices+1,2*NumVertices]];

If[OccSites1 == {} || OccSites2 == {},
       (*cannot associate*)
        Return[{DimerGraph,OccSites,AssocEdges,Valency,Energy}];
  ];

(*the remaining occupied sites on each RNA that can associate*)

OccSites1=Complement[OccSites1,Intersection[OccSites1,VertAssocE]];

OccSites2=Complement[OccSites2,Intersection[OccSites2,VertAssocE]];

randOcc1=OccSites1[[RandomInteger[{1,Length[OccSites1]}]]];

randOcc2=OccSites2[[RandomInteger[{1,Length[OccSites2]}]]];


If[Length[Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc1,1],{randOcc1}],OccSites1]] < Valency && Length[Intersection[Complement[NeighborhoodVertices[DimerGraph,randOcc2,1],{randOcc2}],OccSites2]] < Valency,
     DimerGraph=EdgeAdd[DimerGraph,UndirectedEdge[randOcc1,randOcc2]];
     AssocEdges=Append[AssocEdges,UndirectedEdge[randOcc1,randOcc2]];
     Energy=Energy+epsilon; 
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


StaggerGraph[dimergraph_]:= (

DimerGraph=dimergraph;

Elist=EdgeList[DimerGraph];

NumVertices=Length[VertexList[DimerGraph]];

NewElist={};


(*for every edge, get vertices, stagger them both by the number of vertices in the dimergraph*)
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

(*EdgeUppH[dimergraph_,occsites_,edel_,numedel_,epsilon_,energy_]:= (

DimerGraph=dimergraph;

OccSites=occsites;

Energy=energy;

(*list of deleted edges*)
Edel=edel;

NumEdel=numedel;

(*the removeable edges*)
ERemove=EdgeList[DimerGraph];
sucess=0;
Edel1={};
While[ sucess < NumEdel,
posns={};
	While[Length[posns] < NumEdel,
   	 posns = Union[posns,DeleteDuplicates[RandomSample[Range[Length[ERemove]], NumEdel]]];
    	 ];
(*Now go, edge by edge, to see which connect occupied sites*)
	Do[
	VrandEdge=VertexList[Graph[ERemove[[i]]]];
              (*if both vertices are occupied sites, you gotta use Metropolis to break *)
		If[Length[Intersection[occsites,VrandEdge]] == 2,
                       (*the edge connects two occupied sites, so break with Metro*)
                       (*DelE=0-(epsilon), so -DelE=epsilon*)
                       r= RandomReal[];
                       If[r<Exp[epsilon],
                             (*accept move*)
                              sucess=sucess+1;
                              Edel1=Append[Edel1,ERemove[[i]]];
                              DimerGraph=EdgeDelete[DimerGraph,ERemove[[i]]];
                                  
 (*STOPPED HERE!!!!!!!!!!!!!*)
		  ];
	,{i,posns}];



  ];

posns={};

(*keep track of which edges you will delete*)
Edel1={};

Do[
Edel1=Append[Edel1,ERemove[[i]]];
,{i,posns}];

DimerGraph=EdgeDelete[DimerGraph,Edel1];




) *)

RaisepH[dimergraph_,occsites_,vdel_,edel_,numvdel_]:= (

DimerGraph=dimergraph;

OccSites=occsites;

Vdel=vdel;

Edel=edel;

NumVdel=numvdel;

(*we will randomly remove vertices, but not the ones which are occupied*)

RemoveableVertices=Complement[VertexList[DimerGraph],OccSites];

(*randomly choose numvdel to delete*)
Vdel1={};
While[Length[Vdel1] < numvdel,
   Vdel1 = Union[Vdel1,DeleteDuplicates[RandomSample[Range[Length[RemoveableVertices]],numvdel]]];
      ];

(*for each vertex we delete, we need to keep track of the edges*)
Edel1={};

Do[
(*find neighborhood vertices*)
nhvert=Complement[NeighborhoodVertices[DimerGraph,i,1],{i}];
nhedges={};
	Do[
	nhedges=Append[nhedges,UndirectedEdge[i,j]];
	,{j,nhvert}];

Edel1=Append[Edel1,nhedges];

,{i,Vdel1}];

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

UnOccSites=Complement[Range[Length[VertexList[DimerGraph]]],OccSites];


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


MCMoveGraph[dimergraph_,occsites_,nummoves_,epsilon_,energy_] := (

DimerGraph=dimergraph;

OccSites=occsites;

Energy=energy;

UnOccSites=Complement[Range[Length[VertexList[DimerGraph]]],OccSites];


If[Length[OccSites] == 0 || Length[UnOccSites] == 0,
     Return[{DimerGraph,OccSites,Energy}];
  ];

Do[
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


,{i,1,nummoves}];

Return[{DimerGraph,OccSites,Energy}];

)


MCJumpGraph[dimergraph1_,occsites1_,dimergraph2_,occsites2_,epsilon_,energy_]:=(

DimerGraph1=dimergraph1;

OccSites1=occsites1;

Energy=energy;

UnOccSites1=Complement[Range[Length[VertexList[DimerGraph1]]],OccSites1];

DimerGraph2=dimergraph2;

OccSites2=occsites2;

UnOccSites2=Complement[Range[Length[VertexList[DimerGraph2]]],OccSites2];


(*first we'll choose a random OCCUPIED site on dimergraph1 and an UNOCCUPIED site on dimergraph2*)

If[Length[OccSites1] == 0 || Length[UnOccSites2] == 0,
      Return[{DimerGraph1,OccSites1,DimerGraph2,OccSites2,DimerGraph2,Energy}];
  ];

(*Pick a random occupied site dimergraph1 and a random unoccupied site from dimergraph2*)
     randOcc1=OccSites1[[RandomInteger[{1,Length[OccSites1]}]]];
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

       (*Pick a random OCCUPIED site from DimerGraph2 and a random UNOCCUPIED site from DimerGraph1*)
	If[Length[OccSites2] == 0 || Length[UnOccSites1] == 0,
      	Return[{DimerGraph1,OccSites1,DimerGraph2,OccSites2,DimerGraph2,Energy}];
 	 ];

       randOcc2 = OccSites2[[RandomInteger[{1,Length[OccSites2]}]]];
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


Return[{DimerGraph1,OccSites1,DimerGraph2,OccSites2,DimerGraph2,Energy}];


)


End[ ]

EndPackage[ ]
