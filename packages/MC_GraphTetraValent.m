BeginPackage["MonteCarlo`"]

Needs["GraphUtilities`"];

MCGraphTetraValent::usage = "blah blah"

Begin["`Private`"]

MCGraphTetraValent[adjmat_,numocc_, linkedlist_,numsteps_,numbreaths_,epsilon_] := (
(*Make a graph from the adjacency matrix*)
AdjMatList = adjmat;

AdjGraph = AdjacencyGraph[AdjMatList];
(*number of vertices*)
NumSites = Length[VertexList[AdjGraph]];
(*number of occupied sites*)
NumOcc = numocc;
(*randomly populate the graph*)
OccupiedSites = {};
While[Length[OccupiedSites] < NumOcc,
	OccupiedSites = Union[OccupiedSites,DeleteDuplicates[RandomSample[Range[NumSites],NumOcc]]];
     ];

(*We'll make an adjacecy matrix, but only for the proteins*)


ProtAdjMat=ConstantArray[0,{NumSites,NumSites}];



(*initialize before MC steps*)

Do[

(*First check if the dimer has sites available to bind  4 sites per dimer*)
(*Take the dot product of the row corresponding to the dimer, which is just the number of dimers currently bound to it in our Protein Graph.  If this number is LESS than 4, then it can bind another dimer*)



If[Dot[Part[ProtAdjMat,OccupiedSites[[i]]],Part[ProtAdjMat,OccupiedSites[[i]]]] < 4,
  	(*do something*)
        (*get the occupied nearest-neighbors*)
	OccNN=Intersection[NeighborhoodVertices[AdjGraph,OccupiedSites[[i]],1],Complement[OccupiedSites,{OccupiedSites[[i]]}]];
	
	If[ OccNN == {},
                  (*do nothing*)
                  ProtAdjMat=ProtAdjMat;     ,
		j=1;
		While[(Dot[Part[ProtAdjMat,OccupiedSites[[i]]],Part[ProtAdjMat,OccupiedSites[[i]]]] < 4 ) && (j < Length[OccNN]),
	     	 (*continue going through the while loop until either the dimer has 4 nearest-neighbors or you've come to the end of
                   the occupied nearest-nieghbor list*)

		 If[Dot[Part[ProtAdjMat,OccNN[[j]]],Part[ProtAdjMat,OccNN[[j]]]]< 4,
			(*If the nearest-neighbor has room, you can add*)
			ProtAdjMat[[OccupiedSites[[i]],OccNN[[j]]]]=1;
                        ProtAdjMat[[OccNN[[j]],OccupiedSites[[i]]]]=1;
			j=j+1;
			,
                j=j+1;
                ];
               Print[j]; ];





	     
	           
                   ,

   (*do nothing*)
    ProtAdjMat = ProtAdjMat;		]; ];
Print[i];
,{i,1,NumOcc}];

Join[{ProtAdjMat},{OccupiedSites}]

)

End[ ]

EndPackage[ ]
