BeginPackage["RNAAsGraph`"]

DuplexDistance::usage = "This function is an adaptation of the Li Tai's/Aron's ladder distance function for duplexes.  That is given two duplexes I J, it calculates how many duplexes are between I and J.  To use this function, first convert the *.ct to .bp using ct2bp.sh and bp2duplex.sh to truncate the .bp file into base pairs (i j and j i) which belong to distinct duplexes.  Import this file as a 2 x N matrix in the wrapper m file.  The matrix is called the duplex connection table."

RgGraph::usage = "Returns Rg for a tree graph based on Kramer's formula.  Input is a tree graph."


Begin["`Private`"]

DuplexDistance[duplexpos1_,duplexpos2_,duplexconnectiontable_]:=(

truncatedct=duplexconnectiontable[[All,duplexpos1;;duplexpos2]];
(*number of the unique base pairs*)
numuniquebp=Length[DeleteDuplicates[Flatten[truncatedct]]];

(*number of base pairs*)
numbp=Length[Flatten[truncatedct]];

(*number of duplexes*)
duplexdistance=numuniquebp-0.5*numbp;
	If[MemberQ[truncatedct[[1]],truncatedct[[2,1]]],
          duplexdistance=duplexdistance+1;
	   ];

Return[duplexdistance];
)

RgGraph[graph_] := (

(*make a graph from adjacency matrix*)
(*graph=KirchhoffGraph[kmat];*)
(*number of vertices*)
numverts=Length[VertexList[graph]];
(*get the edge list for the graph*)
listofedges=EdgeList[graph];

product=0;

Do[
(*break on edge at a time*)
brokengraph=Graph[Complement[EdgeList[graph], {listofedges[[i]]}]];
(*Get list of connected components*)
conncomps=ConnectedComponents[brokengraph];
If[Length[conncomps] < 2,
         product = product + Length[conncomps[[1]]];
  ,
  product = product + Length[conncomps[[1]]]*Length[conncomps[[2]]];
  ];

,{i,1,Length[listofedges]}];

rg=N[Sqrt[product]/numverts];
ClearAll[product,brokengraph,listofedges,conncomps];

Return[rg];

)

End[ ]

EndPackage[ ]
