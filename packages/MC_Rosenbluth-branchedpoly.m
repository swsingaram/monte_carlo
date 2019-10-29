BeginPackage["MonteCarloRosenbluthBranchedPolymer`"]

Needs["GraphUtilities`"];

GrabSites::usage = "This program selects an occupied  or an unoccupied sites given by f_i*a or (1-f_i) if the previous site is occupied. Here a is the exp(-e) factor.  Importantly,
                    this function takes as an input the list of weights for all possible outcomes.  In the list of weights, if (i,j is the current coordinate, then the position (i,j+1) is the first element, (i,j-1) is the second, (i+1,j) is the third,
                    and (i-1,j is the fourth. The list of weights has eight elements in total.  The first four and the last four correspond to the above moves for occupied and 
                    unoccupied sites, repectively.  The selected position is returned from the routine ConfigSelect."


Grow::usage = "Here we grow the polymer! We follow the Rosenbluth scheme, adapted for branched polymers with bound particles so that each site is occupied with probability p=Nocc/Nsites.  See Daan Frenkel's 'Understanding Molecular Simulations' chapter 13.2  on configurational biased monte carlo (algorithm 24). This function adds one polymeric unit at a time and returns the updated configuration, and the new rosenbluth factor as well.  In addition, we divide rosenbluth factor by a factor of four for computational convenience--large numbers  take up more memory.  PolyCoors is a list containg the polymer's vital information an element of the list takes the form: {coordinates,occupancy state, vertex index}.  The coordinates are the spatial location on a 2D lattice, the occupancy state is 1 for occupied and 0 for unoccupied, and the vertex index is the index of the vertex in the graph which you must supply from Mathematica."

ConfigSelect::usage = "We select a trial position with prob w(i)/Sum(w(i)), this algorithm is 'Algorithm 41' from Daan Frenkel's book in Appendix J."

ListOfWeights::usage = "Based on the configuration of the polymer, the state of the last vertex added (occupied or not), and the fraction of occupied sites remaining (which you get from GrabSites) we  compute the list of weights for each outcome."

VertexNN::usage = "This function returns the vertices the vertex is connected to on the partially grown polymer."

ListOfWeightsCompExp::usage = "Returns a list of weights from the RAWRBFactor_x_x.txt file.  Remember to replace *^ with e in the file BEFORE using this function.  I usually take care of this by modifying the RAWRBFactor file from the beginning.  We use this function in the competition part of the simulation."

IsOccupied::usage = "Returns whether the site in question is occupied by the polymer.  If it returns {}, site is not occupied by the polymer, if it returns 1 it is and the site is occupied with a protein, if it returns 0 then it is not occupied with a protein."

OccUnOccPosns::usage = "Returns positions of occupied and unoccupied sites from PolyCoors."

Move::usage = "Returns move based on the index i (1-8) which samples eight possible moves (2*4 moves because each of the two states can be in four possible positions): {{0,-1},{0,1},{-1,0},{1,0},{0,-1},{0,1},{-1,0},{1,0}};"

Frontier::usage = "Returns the vertices which you still need to add based on the graph, the partially grown graph, and the index of the last vertex added. INPUT-graph, polymercoordinates, vertexnum"

Rg::usage = "Returns Rg from PolyCoors"

Jump::usage = "DEFECTIVE FUNCTION--WE EXCHANGE PARTICLES USING THE COMPETITION M-FILE.  DO NOT USE THIS FUNCTION!!!!!!!!!!!!!Function for protein transfer betweeen two RNAs dressed with protein.  This function will randomly pick an occupied site on one RNA and an unoccupied site on the other, then determine whether the transfer occurs following the Metropolis rule.  The canonical metropolis rule is based on energy differences and will lead to a bias in this simulation since our states are sampled according to their Rosenbluth weights.  We must accept a move based on the respective rosenbluth factors and their respective binomial coefficients (the binomial coefficient is taken care by selecting for occupied sites and unoccupied sites).  See section 13.2.2 of the Frenkely book for a discussion and use your noggin! Finally this function returns the updated polymer coordinates (which includes which sites are occupied with protein)"

Begin["`Private`"]

Frontier[graph_,polymercoordinates_,vertexnum_]:= (
ClearAll[neighborstoadd];
(*get initial conditions*)
(*find this vertex's position*)
(*best way is to find the pattern match {_,_,vertnum} in level 2 not in the ones that may be in the coordinates*)
VertNumPos=Position[polymercoordinates,{_,_,vertexnum},2][[1,1]];
(*get a list of the neighbors from the completed graph*)
finalneighbors=Complement[NeighborhoodVertices[graph,vertexnum,1],{vertexnum}];
(*get a list of the neighbors from a partially grown graph*)
partialneighbors=VertexNN[polymercoordinates,VertNumPos];
If[ Length[Intersection[partialneighbors,finalneighbors]] < Length[finalneighbors],
	(*return which vertices need to be added, we will have to sort this list
        to ultimately decide which site to add*)
	neighborstoadd=Complement[finalneighbors,partialneighbors];
	ClearAll[finalneighbors,partialneighbors];
	Return[neighborstoadd];
      ,
	(*all neighbors have been added*)
	ClearAll[finalneighbors,partialneighbors,neighborstoadd];
    	Return[{}];
  ];
)

Jump[polyconfigurationRNA1_,polyconfigurationRNA2_,wtinitial_,pathtorbf1_,pathtorbf2_]:= (

ClearAll[PolyCoorsRNA1,PolyCoorsRNA2];

PolyCoorsRNA1=polyconfigurationRNA1;
PolyCoorsRNA2=polyconfigurationRNA2;



(*get occupied and unoccupied sites from both RNAs*)
(*start with RNA 1*)
OccUnOccListRNA1=OccUnOccPosns[PolyCoorsRNA1];
(*list of positions for the occupied sites on RNA1*)
OccListRNA1=OccUnOccListRNA1[[1]];
(*number of occupied sites on RNA1*)
numOccRNA1=Length[OccListRNA1];
(*list of positions for the unoccupied sites on RNA1*)
UnOccListRNA1=OccUnOccListRNA1[[2]];
(*number of unoccupied sites on RNA1*)
numUnOccRNA1=Length[UnOccListRNA1];
(*number of sites on RNA1*)
numbersitesRNA1=numOccRNA1+numUnOccRNA1;

(*now checkout RNA 2*)
OccUnOccListRNA2=OccUnOccPosns[PolyCoorsRNA2];
(*list of positions for the occupied sites on RNA2*)
OccListRNA2=OccUnOccListRNA2[[1]];
(*number of occupied sites on RNA2*)
numOccRNA2=Length[OccListRNA2];
(*list of positions for the unoccupied sites on RNA2*)
UnOccListRNA2=OccUnOccListRNA2[[2]];
(*number of unoccupied sites on RNA2*)
numUnOccRNA2=Length[UnOccListRNA2];
(*number of sites on RNA2*)
numbersitesRNA2=numOccRNA2+numUnOccRNA2;

(**********************************
************************************
REMEMBER DETAILED-BALANCE
***********************************
**********************************
MAYBE JUST FLIP A COIN TO SEE WHO 
WILL BE THE PROTEIN-DONOR
*********************************
********************************
********************************)
l1=RandomReal[];

focc=N[numOccRNA1/(numOccRNA1+numOccRNA2)];

l2=RandomReal[];

funocc=N[numUnOccRNA1/(numUnOccRNA1+numUnOccRNA2)];

(*pick a random occupied site and a random unoccupied site, if they
happen to be on different molecules attempt to transfer*)


If[ l1 <= focc && l2 > funocc,
      (*attempt to transfer from RNA1 to RNA2*)
	(*now we have the positions*)
	If[OccListRNA1 != {} && UnOccListRNA2 != {},
       		 (*so RNA1 has occupied sites and RNA2 has unoccupied sites*)
        	 (*if the transfer is accepted it will put one more particle on RNA2 and remove a particle from RNA1*)
		  (*if the transfer is accepted the RNA1 and RNA2 will have the following particle numbers*)
		numOccRNA1Final=numOccRNA1-1;
		numOccRNA2Final=numOccRNA2+1;
		(*path to file for RNA1*)
		rbfpathRNA1Final=StringJoin[pathtorbf1,"RAWRBFactor_",ToString[numOccRNA1Final],"_",ToString[numbersitesRNA1],".txt"];

		(*path to files for RNA2*)
		rbfpathRNA2Final=StringJoin[pathtorbf2,"RAWRBFactor_",ToString[numOccRNA2Final],"_",ToString[numbersitesRNA2],".txt"];

		
		
		(*pick a configuration for RNA1 at random using the Boltzmann distribution*)
		LOWRNA1=ListOfWeightsCompExp[rbfpathRNA1Final];
		posn=ConfigSelect[LOWRNA1];
		rbfrna1final=LOWRNA1[[posn]];
		wtrna1final=rbfrna1final;
	
		(*pick a configuration for RNA2 at random using Boltzmann distribution*)
		LOWRNA2=ListOfWeightsCompExp[rbfpathRNA2Final];
		posn=ConfigSelect[LOWRNA2];
		rbfrna2final=LOWRNA2[[posn]];
		wtrna2final=rbfrna2final;

        	(*total weight for the configuration RNA1&RNA2*)
		wtfinal=wtrna1final*wtrna2final;
                
		(*if we pick a random real number between 0 and 1 and its less than the minimum of 1 and wtfinal/wtinitial, we accept*)
		r=RandomReal[];
		accratio=wtfinal/wtinitial;
		Print[{wtinitial,wtfinal,accratio,"RNA1 to RNA2"}];
		If[ r < Min[1,accratio],
			(*accept!!*)
			(*new particle distributions*)
			numOccRNA1=numOccRNA1Final;
			numOccRNA2=numOccRNA2Final;	
		 ];	
		

		

(*		(*pick a random position of an occupied site on RNA1*)
		randOccPosn=OccListRNA1[[RandomInteger[{1,Length[OccListRNA1]}]]];
        	(*pick a random position of an unoccupied site on RNA2*)
        	randUnOccPosn=UnOccListRNA2[[RandomInteger[{1,Length[UnOccListRNA2]}]]];
		
		(*calculate the initial energy*)
       		 (*we just need to calculate the number of occupied nearest-neighbors*)
       		 (*get coordinate of the occupied site*)
       		 OccCoor=PolyCoorsRNA1[[randOccPosn,1]];
        	 numocc=0;
       			 Do[
	 			 TempCoor=OccCoor;
         			 TempCoor[[2]]=TempCoor[[2]]+j;
         			 OccupancyState=IsOccupied[PolyCoorsRNA1,TempCoor];
	 			 If[OccupancyState == 1,
                			numocc=numocc+1;
            			   ];
			 ,{j,{-1,1}}];

			Do[
	 			 TempCoor=OccCoor;
          			 TempCoor[[1]]=TempCoor[[1]]+i;
	  			 OccupancyState=IsOccupied[PolyCoorsRNA1,TempCoor];
	  			 If[OccupancyState == 1,
              				 numocc=numocc+1;
           			  ];
			,{i,{-1,1}}];
		 (*calculate initial energy*)
		InitialEnergy = numocc*epsilon;

		(*calculate the final energy*)
		(*get coordinate of the unoccupied site on RNA2*)
        	UnOccCoor=PolyCoorsRNA2[[randUnOccPosn,1]];
		numocc=0;
		Do[
			TempCoor=UnOccCoor;
			TempCoor[[2]]=TempCoor[[2]]+j;
			OccupancyState=IsOccupied[PolyCoorsRNA2,TempCoor];
			If[OccupancyState == 1,
              			numocc=numocc+1;
          	   	  ];
		,{j,{-1,1}}];

		Do[
			TempCoor=UnOccCoor;
			TempCoor[[1]]=TempCoor[[1]]+i;
			OccupancyState=IsOccupied[PolyCoorsRNA2,TempCoor];
        		If[OccupancyState == 1,
               			numocc=numocc+1;
         		  ];
          	,{i,{-1,1}}];
	(*calculate final energy*)
	FinalEnergy = numocc*epsilon;

	DelEnergy = FinalEnergy - InitialEnergy;
		If[DelEnergy <= 0,
                  (*accept move*)
	          (*change occupancy state of RNA1 and RNA2*)
		  PolyCoorsRNA1[[randOccPosn,2]]=0;
		  PolyCoorsRNA2[[randUnOccPosn,2]]=1;
		  Energy = Energy + DelEnergy;
		  ,
                 (*DelEnergy > 0, accept with boltzmann probability*)
                  r = RandomReal[];
			If[r<Exp[-(DelEnergy)],
			   (*Accept Move*)
			   PolyCoorsRNA1[[randOccPosn,2]]=0;
                           PolyCoorsRNA2[[randUnOccPosn,2]]=1;		
			   Energy = Energy + DelEnergy;

			  
                           (*reject move*)
			   
                          ];

		 ];

  	this comments everything except for the first if statement in the true section!*) ];
,

     If[l1 > focc && l2 <= funocc,
          (*attempt transfer from RNA2 to RNA1*)
	If[OccListRNA2 !={} && UnOccListRNA1 != {},
        	(*so RNA2 has occupied sites and RNA1 has unoccupied sites*)	

        	(*if the transfer is accepted it will put one more particle on RNA1 and remove a particle from RNA2*)
		  (*if the transfer is accepted the RNA1 and RNA2 will have the following particle numbers*)
		numOccRNA1Final=numOccRNA1+1;
		numOccRNA2Final=numOccRNA2-1;
		(*path to file for RNA1*)
		rbfpathRNA1Final=StringJoin[pathtorbf1,"RAWRBFactor_",ToString[numOccRNA1Final],"_",ToString[numbersitesRNA1],".txt"];

		(*path to files for RNA2*)
		rbfpathRNA2Final=StringJoin[pathtorbf2,"RAWRBFactor_",ToString[numOccRNA2Final],"_",ToString[numbersitesRNA2],".txt"];
				
		
		(*pick a configuration for RNA1 at random using the Boltzmann distribution*)
		LOWRNA1=ListOfWeightsCompExp[rbfpathRNA1Final];
		posn=ConfigSelect[LOWRNA1];
		rbfrna1final=LOWRNA1[[posn]];
		wtrna1final=rbfrna1final;
	
		(*pick a configuration for RNA2 at random using Boltzmann distribution*)
		LOWRNA2=ListOfWeightsCompExp[rbfpathRNA2Final];
		posn=ConfigSelect[LOWRNA2];
		rbfrna2final=LOWRNA2[[posn]];
		wtrna2final=rbfrna2final;

        	(*total weight for the configuration RNA1&RNA2*)
		wtfinal=wtrna1final*wtrna2final;

		(*if we pick a random real number between 0 and 1 and its less than the minimum of 1 and wtfinal/wtinitial, we accept*)
		r=RandomReal[];
		accratio=wtfinal/wtinitial;
		Print[{wtinitial,wtfinal,accratio,"RNA2 to RNA1"}];
		If[ r < Min[1,accratio],
			(*accept!!*)
			(*new particle distributions*)
			numOccRNA1=numOccRNA1Final;
			numOccRNA2=numOccRNA2Final;	
		 ];	

	(*	(*pick a random position of an occupied site on RNA2*)
        	randOccPosn=OccListRNA2[[RandomInteger[{1,Length[OccListRNA2]}]]];
        	
		(*pick a random position of an unoccupied site on RNA1*)
        	randUnOccPosn=UnOccListRNA1[[RandomInteger[{1,Length[UnOccListRNA1]}]]];
		(*calculate the initial energy*)
        	(*we just need to calculate the number of occupied nearest-neighbors*)
       		 (*get coordinate of the occupied site*)
        	OccCoor=PolyCoorsRNA2[[randOccPosn,1]];
        	numocc=0;
       		 Do[
	  		TempCoor=OccCoor;
          		TempCoor[[2]]=TempCoor[[2]]+j;
          		OccupancyState=IsOccupied[PolyCoorsRNA2,TempCoor];
	  		If[OccupancyState == 1,
                		numocc=numocc+1;
           		 ];
	 	,{j,{-1,1}}];

		Do[
	  		TempCoor=OccCoor;
          		TempCoor[[1]]=TempCoor[[1]]+i;
	  		OccupancyState=IsOccupied[PolyCoorsRNA2,TempCoor];
	  		If[OccupancyState == 1,
               			numocc=numocc+1;
            		 ];
		,{i,{-1,1}}];
		
	
		
		 (*calculate initial energy*)
		 InitialEnergy = numocc*epsilon;
       		(*calculate the final energy*)
		(*get coordinate of the unoccupied site on RNA1*)
        	UnOccCoor=PolyCoorsRNA1[[randUnOccPosn,1]];
		numocc=0;
		Do[
			TempCoor=UnOccCoor;
			TempCoor[[2]]=TempCoor[[2]]+j;
			OccupancyState=IsOccupied[PolyCoorsRNA1,TempCoor];
			If[OccupancyState == 1,
              			numocc=numocc+1;
         		 ];
		,{j,{-1,1}}];

		Do[
			TempCoor=UnOccCoor;
			TempCoor[[1]]=TempCoor[[1]]+i;
			OccupancyState=IsOccupied[PolyCoorsRNA1,TempCoor];
        		If[OccupancyState == 1,
               			numocc=numocc+1;
          		  ];
          	,{i,{-1,1}}];
		(*calculate final energy*)
		FinalEnergy = numocc*epsilon;
		DelEnergy = FinalEnergy - InitialEnergy;
		If[DelEnergy <= 0,
                  	(*accept move*)
	          	(*change occupancy state of RNA1 and RNA2*)
		 	 PolyCoorsRNA2[[randOccPosn,2]]=0;
		  	PolyCoorsRNA1[[randUnOccPosn,2]]=1;
		  	Energy = Energy + DelEnergy;
		  	,
                 	(*DelEnergy > 0, accept with boltzmann probability*)
                  	r = RandomReal[];
				If[r<Exp[-(DelEnergy)],
			   		(*Accept Move*)
			   		PolyCoorsRNA2[[randOccPosn,2]]=0;
                           		PolyCoorsRNA1[[randUnOccPosn,2]]=1;		
			   		Energy = Energy + DelEnergy;

			   	
                           		(*reject move*)
			   		
                          	];

		 ];
     	 this comments everything in the else section *)];

      ];

];
ClearAll[OccUnOccListRNA1,OccListRNA1,UnOccListRNA1,OccUnOccListRNA2,OccListRNA2,UnOccListRNA2,LOWRNA1,LOWRNA2,wtrna2final,wtrna1final,wtfinal];

Return[{numOccRNA1,numOccRNA2}];

)


Rg[polyconfiguration_]:= (


PolyConfig=polyconfiguration;
ChainLength=Length[PolyConfig];
PolyCoors={};

Do[
PolyCoors=Append[PolyCoors,PolyConfig[[j]][[1]]];
,{j,1,Length[PolyConfig]}];

CoorPairs=Subsets[PolyCoors,{2}];
NumCoorPairs=Length[CoorPairs];
RadGy=0;
Do[
RadGy=RadGy+SquaredEuclideanDistance[CoorPairs[[i]][[1]],CoorPairs[[i]][[2]]];
,{i,1,NumCoorPairs}];

ClearAll[PolyConfig,PolyCoors,CoorPairs];
Return[N[Sqrt[RadGy]/ChainLength]];

)

ListOfWeights[lastvertcoordinate_,polycoordinates_,f_,epsilon_]:= (

ClearAll[LOW];

(*let's start sniffing around*)
(*first let's put 0s in the list corresponding to moves that would cause the polymer to run into itself*) 
cycle=0;
FreePosns={};
(*list of weights, constant array of 0s*)
LOW=ConstantArray[0,8];
Do[
TempCoor=lastvertcoordinate;
TempCoor[[2]]=TempCoor[[2]]+j;
cycle=cycle+1;

OccupancyState=IsOccupied[polycoordinates,TempCoor];
If[OccupancyState == {},
    (*this position in LOW is free*)
    FreePosns=Append[FreePosns,cycle];
    (*could ultimately be an unoccupied site, so the corresponding move would be*)      
    FreePosns=Append[FreePosns,cycle+4];
      ,
    (*site is occupied*)
    LOW[[cycle]]=0;
    LOW[[cycle+4]]=0;
    ];
,{j,{-1,1}}];


Do[
TempCoor=lastvertcoordinate;
TempCoor[[1]]=TempCoor[[1]]+i;
cycle=cycle+1;
(*moving up and down *)
OccupancyState=IsOccupied[polycoordinates,TempCoor];
If[OccupancyState == {},
    (*this position in LOW is free*)
    FreePosns=Append[FreePosns,cycle];
    (*could ultimately be an unoccupied site, so the corresponding move would be*)
    FreePosns=Append[FreePosns,cycle+4];
    ,
    (*site is occupied*)
    LOW[[cycle]]=0;
    LOW[[cycle+4]]=0;
  ];
,{i,{-1,1}}];

If[FreePosns=={},
    Return[LOW];
  ];

(*at this point we know where the polymer can move. The list of those positions are in FreePosns*)
(*we will update the list of weights*)




Do[
(*define a temporary coordinate*)
TempCoor=lastvertcoordinate;

If[Posn < 5,
 (*these positions correspond to adding an occupied site*)

(*this corrseponds to the position in ListOfWeight, its the position indexed by Move*)
TempMove=Move[Posn];
TempCoor[[1]]=TempCoor[[1]]+TempMove[[1]];
TempCoor[[2]]=TempCoor[[2]]+TempMove[[2]];

(*now sniff around at new coordinate and keep track of number of occupied sites*)
numocc=0;

	Do[
	TempCoor2=TempCoor;
        TempCoor2[[2]]=TempCoor2[[2]]+j;
	OccupancyState=IsOccupied[polycoordinates,TempCoor2];
	If[OccupancyState == 1,
           numocc=numocc+1;
          ];	

	,{j,{-1,1}}];

	Do[
	TempCoor2=TempCoor;
        TempCoor2[[1]]=TempCoor2[[1]]+i;
	OccupancyState=IsOccupied[polycoordinates,TempCoor2];

        If[OccupancyState == 1,
           numocc=numocc+1;
          ];

         ,{i,{-1,1}}];

Energy=epsilon*numocc;
A=Exp[-Energy];
(*now we need to know how the fraction of occupied sites in the container*)
Weight=A*f;
(*update list of weights*)
LOW[[Posn]]=Weight;
,

(*these positions correspond to adding an unoccupied site*)

LOW[[Posn]]=(1-f);

];

,{Posn,FreePosns}];

ClearAll[FreePosns,TempCoor,OccupancyState,Energy,A,Weight,TempCoor2];
        
Return[LOW];

)

VertexNN[polycoordinates_,posn_]:= (

ClearAll[neighbors];

(*get coordinate corresponding to the position*)
InitialCoor=polycoordinates[[posn,1]];
(*now start sniffing around using the code adapted from list of weights*)
(*now sniff around at new coordinate and keep track of the sites occupied by the polymer*)
neighbors={};

	Do[
	TempCoor=InitialCoor;
        TempCoor[[2]]=TempCoor[[2]]+j;
	OccupancyState=IsOccupied[polycoordinates,TempCoor];
	If[NumberQ[OccupancyState],
	   (*get position of the coordinate in the list polycoordinates*)
	   location=Position[polycoordinates,TempCoor];
	   pos=polycoordinates[[location[[1,1]],3]];
          neighbors=Append[neighbors,pos];
          ];	

	,{j,{-1,1}}];

	Do[
	TempCoor=InitialCoor;
        TempCoor[[1]]=TempCoor[[1]]+i;
	OccupancyState=IsOccupied[polycoordinates,TempCoor];
        If[NumberQ[OccupancyState],
	   (*get position of the coordinate in the list polycoordinates*)
	   location=Position[polycoordinates,TempCoor];
	   pos=polycoordinates[[location[[1,1]],3]];
           neighbors=Append[neighbors,pos];
          ];

         ,{i,{-1,1}}];
ClearAll[InitialCoor,TempCoor,OccupancyState,location,pos];
Return[neighbors];
)

ListOfWeightsCompExp[rawrbfactorfile_]:= (
(*import data from the RBFactor file*)
LOW=Import[rawrbfactorfile,"List"];
Return[LOW];
)



IsOccupied[polycoordinates_,coordinate_]:= (
(*get position of the coordinate in the list polycoordinates*)
pos=Position[polycoordinates,coordinate];
If[pos == {},
     (*site is not occupied by the polymer*)
     OccupancyState=pos;
     ClearAll[pos];
     Return[OccupancyState];
  ];

(*change position second element to 2 to get at occupancy state*)
pos[[1,2]]=2;
OccupancyState=polycoordinates[[pos[[1]][[1]],pos[[1]][[2]]]];
ClearAll[pos];
Return[OccupancyState];
)

OccUnOccPosns[polymercoordinates_] := (
ClearAll[OccPosns,UnOccPosns];

OccPosns={};
UnOccPosns={};
numsites=Length[polymercoordinates];
Do[
OccupancyState=polymercoordinates[[i,2]];
(*Coordinate=polymercoordinates[[i,1]];*)

If[OccupancyState == 1,
        (*site is occupied with protein*)
        (*keep track of its position*)
        OccPosns=Append[OccPosns,i];
        ,
       (*sites is unoccupied*)
       UnOccPosns=Append[UnOccPosns,i];
  ];
,{i,1,numsites}];

Return[{OccPosns,UnOccPosns}];

)


GrabSites[numocc_,numunocc_,position_]:= (
(*numocc and numunocc are the number of occupied sites and unoccupied sites left*)
NumOcc=numocc;
NumUnOcc=numunocc;

If[position < 5,
  (*site is occupied*)
  NumOcc=NumOcc-1;,
  
  (*site is unoccupied*)
  NumUnOcc=NumUnOcc-1;
  ];

Return[{NumOcc,NumUnOcc}];

)

ConfigSelect[listofweights_]:= (
(*the sites sack contain occupied sites and unoccupied sites, represented by 1s and 0s.*)
(*pick a random real number between 0 and 1*)

r=RandomReal[];
(*now implement the select subroutine*)
(*first compute the random cutoff*)
CutOff=N[Total[listofweights]]*r;
CumeZ=listofweights[[1]];
(*keep increasing CumeZ until we are passed the cutoff*)
i=1;
While[CumeZ < CutOff,
i=i+1;
CumeZ=CumeZ+listofweights[[i]];
];
(*return the value i*)
ClearAll[CutOff,CumeZ,r];
Return[i];
(*This is the position.Make sure you update the distribution of occupied and unoccupied sites at some point, use GrabSites*)

)

Move[position_] := (

MovesList={{0,-1},{0,1},{-1,0},{1,0},{0,-1},{0,1},{-1,0},{1,0}};

TempMove=MovesList[[position]];

Return[TempMove];

)



Grow[polycoordinates_,position_,vertexindex_,lastvertcoordinate_,listofweights_,w_,f_]:= (
(*we're assuming you can move, we'll return the updated configuration along with the rosenbluth factor*)
ClearAll[PolyCoors,RBFactor,NextVertCoor,DelE];

PolyCoors=polycoordinates;


(*list of moves based on position*)


LocalZ=N[Total[listofweights]];

(*updated rosenbluth factor, put in the factor of 4 thing*)
RBFactor=(LocalZ*w)/4;
(*determine the next coordinate*)
 AcceptedMove=Move[position];
 NextVertCoor = {lastvertcoordinate[[1]]+AcceptedMove[[1]],lastvertcoordinate[[2]]+AcceptedMove[[2]]};

If[position < 5,
  
(*site is occupied*)
 NextPolyCoors = {NextVertCoor, 1, vertexindex};
 DelE=-Log[N[listofweights[[position]]/f]];
,

 (*site is unoccupied*)
 NextPolyCoors = {NextVertCoor, 0, vertexindex};
 DelE=-Log[N[listofweights[[position]]/(1-f)]];
 ];
PolyCoors=Append[PolyCoors,NextPolyCoors];
(*return the updated rosenbluth factor and the updated polymer configuration*)

ClearAll[NextPolyCoors,LocalZ];

Return[{PolyCoors, NextVertCoor,RBFactor,DelE}];

)



End[ ]

EndPackage[ ]
