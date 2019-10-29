BeginPackage["test`"]

Needs["GraphUtilities`"];

MCBindAndMoveGraph::usage = "sfdasf"

Begin["`Private`"]

MCBindAndMoveGraph[adjmatfile_, occupiedsites_,numsteps_,numbreaths_,epsilon_] := (
(*adjmatfile is the RAG ourput.  It contains Degree Matrix, so make sure you take this out first. Remove the brackets and commas too.*)
AdjMatList = Import[adjmatfile,"Table"];
(*Make Tree graph from adjacency matrix*)
AdjGraph = AdjacencyGraph[AdjMatList];
(*check to see if tree graph is connected*)
If[TreeGraphQ[AdjGraph] == False, Return["Tree graph representation of secondary structure is not connected."];];

no = Length[occupiedsites];

(*File names*)
occsitesfilename = StringJoin[ToString[no], "_OccSitesconfig.dat"];

(*vertices of the unoccupied sites*)
UnoccupiedSites = Complement[VertexList[AdjGraph],occupiedsites];

(*Pick a random unoccupied site, this site is now occupied*)
RandomOcc = UnoccupiedSites[[RandomInteger[{1,Length[UnoccupiedSites]}]]];
OccupiedSites = Append[occupiedsites,RandomOcc];
UnoccupiedSites = Complement[UnoccupiedSites,{RandomOcc}];
TaggedOcc = RandomOcc;
(*Metropolois Algorithm*)
Energy=0;
Do[

(*
********************************
********************************
FIRST ALLOW DIMER TO MOVE ALONG
 THE RNA
********************************
********************************
*)

(*figure out the unoccupied nearest-neighbors*)
UnOccNN = Intersection[Complement[NeighborhoodVertices[AdjGraph,TaggedOcc,1],{TaggedOcc}],UnoccupiedSites];
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
];];

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
   (*Using TaggedOcc, figure out the occupied sites TWO links away*)
   OccSitesTwoLinks = Intersection[Complement[NeighborhoodVertices[AdjGraph,TaggedOcc,2],{TaggedOcc}],OccupiedSites];
   RandOcc = OccSitesTwoLinks[[RandomInteger[{1,Length[OccSitesTwoLinks]}]]]; 
(*Determine if we allow the two proteins to interact by treating
the three-vertex subgraph as an ideal polymer and accepting the move
based on end-to-end ideal polymer statistics.*)
(*So I'll accept the move if the end-to-end distance is less than or equal to 2 nm.  I've assumed a duplex is 1.5 nm (5bp and 0.3 nm/bp).*)   
   (*The probability ends up being 0.46*)
r = RandomReal[];
    If[r <= 0.46,
        (*Accept Move *)
        (*Add link between two vertices*) 
         AdjMatList[[TaggedOcc,RandOcc]]=1;
         AdjMatList[[RandOcc,TaggedOcc]]=1;
         AdjGraph=AdjacencyGraph[AdjMatList];
        ,
        
        (*Reject Move*)
        AdjGraph = AdjGraph;
       
          
        
      ]; 
 
  ];


Print[Energy];

,{i,1,numsteps}]

Print[OccupiedSites];

)

End[ ]

EndPackage3
[ ]
