Needs["MonteCarloRosenbluthBranchedPolymer`","/u/home/m/moosee33/scripts/mathematica/packages/MC_Rosenbluth-branchedpoly.m"];
ParallelEvaluate[Needs["MonteCarloRosenbluthBranchedPolymer`","/u/home/m/moosee33/scripts/mathematica/packages/MC_Rosenbluth-branchedpoly.m"]];
(*import the compact tree*)
rnaadjmat=Import["/u/home/m/moosee33/data/lapmats/BMV50tree.lap","Table"];
rnaasgraph=AdjacencyGraph[rnaadjmat];
numsites=Length[VertexList[rnaasgraph]];
epsilon=-3;
$HistoryLength=1;
(*let us try to grow a branched polymer*)
(*number of cp*)
DistributeDefinitions[rnaasgraph,epsilon,numsites,$HistoryLength];
TimeConstrained[
ParallelDo[
(*number of configurations*)
numruns=10;
(*initialize while loop*)
i=1;
If[ FileExistsQ[StringJoin["radiusgyration_",ToString[cp],"_",ToString[numsites],".txt"]],
	(*we will skip cp,rna ratios that we've done already*)	
	numlines=First[ReadList[StringJoin["!cat radiusgyration_",ToString[cp],"_",ToString[numsites],".txt"," | wc -l"]]];
	If[ numlines >= numruns ,
     		i=numruns+1;
  	  ];
   ];
Print[{cp,numsites}];
While[i <= numruns,
	        
	Unprotect[In,Out];
	Clear[In,Out];
	Protect[In,Out];
	ClearSystemCache[];
	ClearAll[r, FrontMatrix,f,NumOcc,NumUnOcc,CPRNAratio,grablist, PolyCoors,RBFactor,NNlist,LOW,GrowList,Col,VertStem,LastVertCoor,PosnState];
	(*ClearAll[GrabSites,Grow,ConfigSelect,ListOfWeights,VertexNN,IsOccupied,Move,Frontier,Rg];
	*)
         

	(*generate a number between 0 and 1*)
	r=RandomReal[];
	Energy=0;
	(*create an adjacency matrix of the vertices on the frontier*)
	FrontMatrix=ConstantArray[0,{numsites,numsites}];
	f=cp/numsites;
	NumOcc=numsites*f;
	NumUnOcc=numsites*(1-f);
	CPRNAratio=N[f];

	(*initialize the vertex*)
	InitialVertex=RandomInteger[{1,numsites}];


	If[r < CPRNAratio,
		(*randomly selected an occupied site*)
		grablist=GrabSites[NumOcc,NumUnOcc,1];
		NumOcc=grablist[[1]];
		NumUnOcc=grablist[[2]];
		(*initialize coordinates*)
		PolyCoors={{{0,0},1,InitialVertex}};
		f=N[NumOcc/(NumOcc+NumUnOcc)];

		,

		(*randomly selected an unoccupied site*)	
		grablist=GrabSites[NumOcc,NumUnOcc,5];
		NumOcc=grablist[[1]];
		NumUnOcc=grablist[[2]];
		(*initialize coordinates*)
		PolyCoors={{{0,0},0,InitialVertex}};
		f=N[NumOcc/(NumOcc+NumUnOcc)];

	];

Result=Catch[DepthFirstScan[rnaasgraph,InitialVertex,{"PrevisitVertex" -> ( 
	VertIndex=#1;
	If[VertIndex == InitialVertex,
	     RBFactor=1;
	     NNlist=Frontier[rnaasgraph,PolyCoors,InitialVertex];
	    Do[
	    (*place ones in front matrix indicating which vertices need to be added to the partially grown graph*)
	    FrontMatrix[[InitialVertex,k]]=1;
	    ,{k,NNlist}];
	   
	  ,
	  (*vertex scanned by depth of scan will definitely be a member of the frontier...i hope*)   
	  NNlist={};
          (*we need to know which vertex, amongst the added vertices, to add to*)
	  (*it is easy to figure this out by searching for the "1" in the column of the frontier matrix*)
	  (*grab the column corresponding to the vertex index*)	    
	  Col=Take[FrontMatrix,All,{VertIndex}];
	  (*vertex stem is the number/position of where the vertex should be added, you can get its coordinates from PolyCoors*)	
          VertStem=Position[Col,{1}][[1,1]];
	  (*get position of vertex stem in polycoors*)
	  VertStemPos=Position[PolyCoors,{_,_,VertStem},2][[1,1]];
	  (*get the vertex coordinate corresponding to the VertStem*)
	  LastVertCoor=PolyCoors[[VertStemPos,1]];
          (*calculate list of weights*)
	  LOW=ListOfWeights[LastVertCoor,PolyCoors,f,epsilon];
	  (*check if LOW has a 0 in every position, this means the polymer got stuck*)
	  If[Total[N[LOW]] == 0,
		(*return empty polycoors and 0 for rosenbluth factor*)
         ClearAll[r, FrontMatrix,f,NumOcc,NumUnOcc,CPRNAratio,grablist, PolyCoors,RBFactor,NNlist,LOW,GrowList,Col,VertStem,LastVertCoor,PosnState]; 
		(*ClearAll[GrabSites,Grow,ConfigSelect,ListOfWeights,VertexNN,IsOccupied,Move,Frontier,Rg];*)
		PolyCoors={};
		RBFactor=0;
		Unprotect[In,Out];
		Clear[In,Out];
		Protect[In,Out];
		ClearSystemCache[];
		Throw[2];
		(*Return[{PolyCoors,RBFactor}];*)
	     ];
	  
	   (*select an outcome: specify posn and occupancy state*)
	   PosnState=ConfigSelect[LOW];

	  (*at this point we can add a new coordinate and occupancy to PolyCoors, as well as update NumOcc,
	   NumUnOcc*)
	  GrowList=Grow[PolyCoors,PosnState,VertIndex,LastVertCoor,LOW,RBFactor,f];

	  (*update polymer coordinates*)
	  PolyCoors=GrowList[[1]];
	(*the vertex has been added, so we will need to update the Frontier matrix*)
	(*FrontMatrix[[VertStem,VertIndex]]=0;*)
	(*now get the frontier vertices for the vertex we just added*)
	(*only update the matrix if there are frontier vertices to add*)
	NNlist=Frontier[rnaasgraph,PolyCoors,VertIndex];
     	If[ Length[NNlist] > 0,
		Do[
		FrontMatrix[[VertIndex,k]]=1;
	 	,{k,NNlist}];
       	   ];
	
	(*update rosenbluth factor*)
	RBFactor=GrowList[[3]];
	(*update energy*)
	DelE=GrowList[[4]];
	Energy=Energy+DelE;
	(*update NumOcc, NumUnOcc*)
	grablist=GrabSites[NumOcc,NumUnOcc,PosnState];
	NumOcc=grablist[[1]];
	NumUnOcc=grablist[[2]];
	(*update fraction of occupied sites*)   
	If[ (NumOcc+NumUnOcc) != 0,
	f=N[NumOcc/(NumOcc+NumUnOcc)];
   	  ];

	  ];

	&)}];];

(*If[Result == 2,
      Print[{cp,numsites}];
  ];*)

If[RBFactor != 0,
PutAppend[PolyCoors,StringJoin["RAWPolyCoors_",ToString[cp],"_",ToString[numsites],".txt"]];
PutAppend[Rg[PolyCoors],StringJoin["radiusgyration_",ToString[cp],"_",ToString[numsites],".txt"]];
PutAppend[RBFactor,StringJoin["RAWRBFactor_",ToString[cp],"_",ToString[numsites],".txt"]];
PutAppend[Energy,StringJoin["RAWEnergy_",ToString[cp],"_",ToString[numsites],".txt"]];
(*increment i*)
i=i+1;
  ];


  ];

,{cp,0,numsites}}];
, 180];

Exit[]
