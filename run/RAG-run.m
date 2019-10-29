Needs["RNAAsGraph`","~/scripts/mathematica/packages/RAG.m"];
(*import the duplex file.  you need to transpose it because
then it will be compatible with the function DuplexDistance*)
pathtofile1="~/test/output_1.bp.duplex"
duplexct=Transpose[Import[pathtofile1,"Table"]];
(*also import the distinct duplexes using their base-pairs
as their index*)
pathtofile2="~/test/output_1.bp.half.duplex"
duplexbps=Transpose[Import[pathtofile2,"Table"]];
(*you only need the first row*)
duplexlist=duplexbps[[1]];
(*all pairs of duplexes, which you will eventually
search for all pairs that are separated by a distance of 2
it is a VERY SLOW way for doing things!*)
numduplexes=Length[duplexlist];
duplexpairs=Subsets[Range[1,numduplexes],{2}];
(*connectionlist=ConstantArray[{},numduplexes];*)
(*In mathematica 10 you could construct the "adjacency matrix" then find the spanning tree
through FindSpanningTree*)
adjmat=ConstantArray[0,{numduplexes,numduplexes}];
testvec=ConstantArray[0,{numduplexes,1}];
(*to make this process a bit faster you can use parallel computing
just make the connectionlist a shared variable*)
Do[
(*get the indices of the two duplexes, neccessarily pair[[1]] < pair[[2]]*)
{pairpos1,pairpos2}={pair[[1]],pair[[2]]};
(*get the first base pair of duplex #1 and duplex #2*)
duplex1=duplexbps[[1,pairpos1]];
duplex2=duplexbps[[1,pairpos2]];
(*now get their positions in the duplex connection table*)
pos1=Position[duplexct,duplex1][[1,2]];
pos2=Position[duplexct,duplex2][[1,2]];
(*now calculate the distance between the two duplexes*)
distance=DuplexDistance[pos1,pos2,duplexct];
(*when DuplexDistance returns 2, it means the duplexes are next to each other*)
If[distance == 2,
   (*connectionlist[[pairpos1]]=Append[connectionlist[[pairpos1]],pairpos2];*)
   (*to update the adjacency matrix, if both vertices (duplexes) have been added, don't update the adjacency matrix*)
   col1=Take[adjmat,All,{pairpos1}];
   col2=Take[adjmat,All,{pairpos2}];
	If[col1 == testvec || col2 == testvec,
			adjmat[[pairpos1,pairpos2]]=1;
			adjmat[[pairpos2,pairpos1]]=1;
	  ];
   
  ];
,{pair,duplexpairs}];
(*Print[connectionlist];*)
Export[StringJoin[pathtofile1,".adjmat"],adjmat,"Table"];
(*at this point we are ready to test!*)
(*connection list is basically the adjacency matrix*)
(*now we'll make the tree, the rule is simply: if both vertices (duplexes) 
have already been placed, you can't add a new connection between them*)
(*calculate Rg*)
rg=RgGraph[AdjacencyGraph[adjmat]];
Print[rg];
