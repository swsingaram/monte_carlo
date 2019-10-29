Needs["MonteCarloRosenbluthBranchedPolymer`","~/scripts/mathematica/packages/MC_Rosenbluth-branchedpoly.m"];
$HistoryLength=1;
(*initialize the number of bound protein on both RNAs*)
numsitesRNA1=50;
numCPRNA1=35;
numsitesRNA2=50;
numCPRNA2=15;
(*numsitesRNA3=26;
numCPRNA3=13;*)

PathToFile1="~/data/compexps/branched/compvsext/compact/3kT/"
PathToFile2="~/data/compexps/branched/compvsext/ext/3kT/"
numsteps=100000;
TimeConstrained[
Do[

	(*ClearAll[ListOfWeightsCompExp,ConfigSelect,Jump,OccUnOccPosns];*)
          

	(*path to files for RNA1*)
	rbfpathRNA1=StringJoin[PathToFile1,"RAWRBFactor_",ToString[numCPRNA1],"_",ToString[numsitesRNA1],".txt"];

	(*path to files for RNA2*)
	rbfpathRNA2=StringJoin[PathToFile2,"RAWRBFactor_",ToString[numCPRNA2],"_",ToString[numsitesRNA2],".txt"];
	
	(*(*path to files for RNA3*)
	rbfpathRNA3=StringJoin[PathToFile3,"RAWRBFactor_",ToString[numCPRNA3],"_",ToString[numsitesRNA3],".txt"];
	polycoorspathRNA3=StringJoin[PathToFile3,"RAWPolyCoors_",ToString[numCPRNA3],"_",ToString[numsitesRNA3],".txt"];*)

	(*pick a configuration for RNA1 at random using the Boltzmann distribution*)
	LOWRNA1=ListOfWeightsCompExp[rbfpathRNA1];
	posn=ConfigSelect[LOWRNA1];
	rbf1=LOWRNA1[[posn]];
	wt1=rbf1*Binomial[numsitesRNA1,numCPRNA1];
	
	(*pick a configuration for RNA2 at random using Boltzmann distribution*)
	LOWRNA2=ListOfWeightsCompExp[rbfpathRNA2];
	posn=ConfigSelect[LOWRNA2];
	rbf2=LOWRNA2[[posn]];
	wt2=rbf2*Binomial[numsitesRNA2,numCPRNA2];

        (*total weight for the configuration RNA1&RNA2*)
	wtinitial=wt1*wt2;
	
	(*now decide whether particle is transferred from RNA1 to RNA2 or
          from RNA2 to RNA1*)

	donor=ConfigSelect[{1,1}];
	(*if donor is equal to 1 then RNA1 transfers to RNA2,
          if donor is equal to 2 then RNA2 transfers to RNA1*)
	If[ donor == 1,
		If[ numCPRNA1 > 0,
		(*RNA1 transfers to RNA2*)
		numCPRNA1Final=numCPRNA1-1;
		numCPRNA2Final=numCPRNA2+1;
		,
		(*can't move*)
		numCPRNA1Final=numCPRNA1;
		numCPRNA2Final=numCPRNA2;
		];
		,
		
		If[ numCPRNA1 < 50,
		(*RNA2 transfers to RNA1*)
		numCPRNA1Final=numCPRNA1+1;
		numCPRNA2Final=numCPRNA2-1;
		,
		(*can't move*)
		numCPRNA1Final=numCPRNA1;
		numCPRNA2Final=numCPRNA2;
		 ];
	];

	(*path to files for RNA1*)
	rbfpathRNA1Final=StringJoin[PathToFile1,"RAWRBFactor_",ToString[numCPRNA1Final],"_",ToString[numsitesRNA1],".txt"];

	(*path to files for RNA2*)
	rbfpathRNA2Final=StringJoin[PathToFile2,"RAWRBFactor_",ToString[numCPRNA2Final],"_",ToString[numsitesRNA2],".txt"];
	
	(*pick a configuration for RNA1 at random using the Boltzmann distribution*)
	LOWRNA1=ListOfWeightsCompExp[rbfpathRNA1Final];
	posn=ConfigSelect[LOWRNA1];
	rbf1final=LOWRNA1[[posn]];
	wt1=rbf1final*Binomial[numsitesRNA1,numCPRNA1Final];
	
	(*pick a configuration for RNA2 at random using Boltzmann distribution*)
	LOWRNA2=ListOfWeightsCompExp[rbfpathRNA2Final];
	posn=ConfigSelect[LOWRNA2];
	rbf2=LOWRNA2[[posn]];
	wt2=rbf2*Binomial[numsitesRNA2,numCPRNA2Final];

        (*total weight for the configuration RNA1&RNA2*)
	wtfinal=wt1*wt2;
	(*accept the move based on their relative weights!*)
	r=RandomReal[];
	accratio=wtfinal/wtinitial;
	If [ r < Min[1,accratio],
		(*accept*)
		(*new particle distributions*)
		numCPRNA1=numCPRNA1Final;
		numCPRNA2=numCPRNA2Final;
	    ]

If[Mod[i,100] == 0,
ClearSystemCache[];
PutAppend[numCPRNA1,"./numcprna1.dat"];
PutAppend[numCPRNA2,"./numcprna2.dat"];
 ];
ClearAll[wtfinal,wtinitial,accratio,LOWRNA1,LOWRNA2,rbf1,rbf2,wt1,wt2];

,{i,1,numsteps}];
,180];












