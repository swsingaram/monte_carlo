(* ::Package:: *)

Off[General::compat]
<< Combinatorica`;
mld[seq_] :=(System`GraphDiameter[System`AdjacencyGraph[Combinatorica`ToAdjacencyMatrix[Combinatorica`CodeToLabeledTree[seq]]]]);

Rg[seq_] := 
  Module[{eli, eln, prodl, tl, ed, new, ll, rr, pr, tprod,ti},
	t= System`AdjacencyGraph[Combinatorica`ToAdjacencyMatrix[Combinatorica`CodeToLabeledTree[seq]]];
   eli = System`EdgeList[t]; eln = Length[eli]; prodl = {};
   tl = VertexCount[t];
   Do[ed = Take[eli, {j}]; new = EdgeDelete[t, ed]; 
    {ll, rr} = System`ConnectedComponents[new];
    pr = Length[ll]*Length[rr]; 
    prodl = Append[prodl, pr]; 
    , {j, eln}
    ];
   tprod = Total[prodl]; 
   N[Round[Sqrt[tprod/(tl^2)],10^-5]]
   
   ];

VertexPosns[plist_] :=
Module[{posns,l},	
posns={};
	l=Length[plist];
	Do[posns=Append[posns,Flatten[Position[plist,i]-l-1]];,{i,Range[Max[plist]]}];
    posns];

(*# find location of high-fold vertices connected to leaves in pruefer sequence*)
LeafLocns[plist_,vposns_] :=
Module[{l,llocns,nnltlocn,nnlt,nnltrtmostloc},
l=Length[plist];
llocns={};
Do[

   nnltlocn=i-1;
   nnlt=plist[[nnltlocn]];
   nnltrtmostloc = Max[vposns[[nnlt]]];

  If[ nnltlocn < nnltrtmostloc,
      llocns=Append[llocns,i];
    ];

,{i,-1,-l+1,-1}];
  llocns=Append[llocns,-1];
  llocns=Append[llocns,-l];
llocns
];

PruferDistanceFast[plist_,vposns_,llocns_]:=
Module[{d,v1,pos,dist},
(*remember to delete -1 from leaf locations if the first position is degree = 2 *)
d=0;
pos=-1;
Catch[While[True,
v1=plist[[pos]];
(*pick a random position in the prufer sequence for this vertex*)
pos=RandomChoice[vposns[[v1]]];
(*did you run into a leaf?*)
If[MemberQ[llocns,pos],
(*end of the line, exit loop*)
d=d;
Throw[d];
];
(*you are next to a deg. 2 vertex or higher, so increase the distance*)
d=d+1;
pos=pos-1;
(*now repeat*)
];
]
];

PruferDistanceLeafLeaf[plist_,vposns_,llocns_]:=
Module[{d,v1,pos,dist},
(*remember to delete -1 from leaf locations if the first position is degree = 2 *)
d=0;
(*choose a random leaf*)
pos=RandomChoice[llocns];
Catch[While[True,
v1=plist[[pos]];
(*move to the highest index of this vertex in the prufer sequence for this vertex*)
pos=Max[vposns[[v1]]];
If[pos==-1,
   (*end of the line, exit loop*)
   d=d;
Throw[d];
   ];
(*move rightwards*)
pos=pos+1;
(*increase distance*)
d=d+1;
(*now repeat*)
];
]
];

PruferDistanceFastBiased[plist_,vposns_,llocns_]:=
Module[{d,v1,pos,posset,dist,possiblepos},
(*remember to delete -1 from leaf locations if the first position is degree = 2 *)
d=0;
pos=-1;
dist=Catch[While[True,
v1=plist[[pos]];
(*pick a random position in the prufer sequence for this vertex*)
posset=Complement[vposns[[v1]],Intersection[llocns,vposns[[v1]]]];
(*did you run into a leaf?*)
If[posset=={},
(*end of the line, exit loop*)
d=d;
Throw[d];
];
pos=RandomChoice[posset];
(*you are next to a deg. 2 vertex or higher, so increase the distance*)
d=d+1;
pos=pos-1;
(*now repeat*)
];
];
dist
];


(*# given two indices i, j (|i| < |j|, prufer sequence, and vertex posns, this
# function finds where the corresponding values s_i and s_j occur so the 
# distance can be calculated.  More explicitly, once we know s_i and s_j we look
# to the right of s_i to find the first occurence of s_j (call it s_j^r)
# and then move left-wards to find the first occurence of s_i (call it s_i^l)*)
FindOccurences[i_,j_,plist_,vposns_]:=
Module[{s1,s2,s1locns,s2locns,s2rlocn,out,s1llocn},
    (*get values of i,j, |i| < |j|*)
    s1=plist[[i]];
    s2=plist[[j]];
    (*locations of s_ 1 and s_ 2*)
    s1locns=vposns[[s1]];
    s2locns=vposns[[s2]];
    (*find position of s_ 2 which is closest to s_ 1, but is still to the right;*)
    s2rlocn ={};
	Do[
	If[ k > i,
        s2rlocn=Append[s2rlocn,k];
      ];
      ,{k,s2locns}];
 
    If[ s2rlocn == {},
        (*there isn't an instance of s_ 2 to the right of s_ 1*)
        out={i,j};,

        s2rlocn = Min[s2rlocn];
        (*find position of s_ 1 which is closest to s_ 2_r, but is to the left;*)
        s1llocn={};
		Do[
		If[ k < s2rlocn,
           s1llocn=Append[s1llocn,k]];
		 ,{k,s1locns}];
		s1llocn=Max[s1llocn];
        out={s1llocn,s2rlocn};
        ];
out
];

(*# now we define the main function of this module --pruefer_distance _list.  Given two
# indices i,j (|i| < |j|), the prufer sequence, and the vertex posns returns
# the vertices involved in the graph distance between them.  
#
# this function is defined recursively*)

PrueferDistanceList[i_,j_,plist_,vposns_]:=
Module[{d,nnleftloc,occurences},
    (*keep track of vertices encountered when moving from i to j*)
    d={};
    Do[
        nnleftloc=k-1;
        occurences=FindOccurences[k, k-1, plist, vposns];
        If[ occurences == {k,k-1},
				If[plist[[k]]==plist[[k-1]],
                      d=d;,
                     d=Append[d,plist[[k-1]]];];,
           d=Append[d,PrueferDistanceList[occurences[[2]],occurences[[1]],plist,vposns]];];
	  ,{k,i,j+1,-1}];    
	d
 ];
(*# now we calculate the distance between two vertices using output from
# pruefer_dist _list  the leaf-leaf distance is 2 more*)  
PrueferDistance[dlist_]:=
Module[{pdistance,distance},
    pdistance=Flatten[dlist];
    distance = Flatten[Position[Tally[pdistance],{_,_Integer?OddQ},2]];
Length[distance]
];

PrueferDistance2[i_,j_,plist_,vposns_,vmet_,d_]:=
Module[{venc,dist,nnleftloc,occurences},
    (*keep track of vertices encountered when moving from i to j*)
  venc=vmet;
  dist=d;
    Do[
        nnleftloc=k-1;
        occurences=FindOccurences[k, k-1, plist, vposns];
        If[ occurences == {k,k-1},
				If[plist[[k]]==plist[[k-1]],
                      dist=dist;,
					If[MemberQ[venc,plist[[k-1]]],
						 venc=Complement[venc,{plist[[k-1]]}];
						 dist=dist-1;,
                         venc=Append[venc,plist[[k-1]]];
						 dist=dist+1;];];,
           {dist,venc}=PrueferDistance2[occurences[[2]],occurences[[1]],plist,vposns,venc,dist]];
	  ,{k,i,j+1,-1}];    
	{dist,venc} ];
(*this matrix constains distances between all neighboring elements in pruefer*)
PruferGraphAuxMatrix[plist_,vposns_]:=
Module[{mat, occurences},
mat=1-IdentityMatrix[Max[plist]];
Do[
occurences=FindOccurences[i,i-1,plist,vposns];
If[occurences!={i,i-1},
mat[[plist[[i]],plist[[i-1]]]]=PruferGraphAuxDistance[occurences[[2]],occurences[[1]],mat,plist,vposns];
];
,{i,-1,-Length[plist]+1,-1}];
mat
];
(*calculate the distance between two vertices*)
PruferGraphAuxDistance[i_,j_,mat_,plist_,vposns_]:=
Module[{d,v1,v2,step,occur},
d=0;
Do[
v1=plist[[k]];
v2=plist[[k-1]];
step=mat[[v1,v2]];
occur=FindOccurences[k,k-1,plist,vposns];
If[occur!={k,k-1},
d=d-step;,
d=d+step];
,{k,i,j+1,-1}];
Abs[d]];

