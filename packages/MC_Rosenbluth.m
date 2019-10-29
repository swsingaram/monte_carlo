BeginPackage["MonteCarloRosenbluth`"]

GrabSites::usage = "This program selects an occupied  or an unoccupied sites given by f_i*a or (1-f_i) if the previous site is occupied. Here a is the exp(-e) factor.  Importantly,
                    this function takes as an input a list of weights.  Let {i,j} be the current coordinate, then {i,j+1} is the first element, {i,j-1} is the second, {i+1,j} is the third,
                    and {i-1,j} is the fourth. Taken as an input are a list of weights which has eight elements.  The first four and the last four correspond to the above moves for occupied and 
                    unoccupied sites, repectively.  The position is given from ConfigSelect."


Grow::usage = "grow a polymer, following the Rosenbluth scheme, following a particular degree distribution such that each site is occupied with probability p=Nocc/Nsites.  See Daan Frenkel's 'Understanding Molecular Simulations' chapter 13.2  on configurational biased monte carlo. This function adds one polymeric unit, returns the updated configuration, and returns rosenbluth factor as well."

ConfigSelect::usage = "We select a trial position with prob w(i)/Sum(w(i)), this algorithm is 'Algorithm 41' from Daan Frenkel's book in Appendix J."

ListOfWeights::usage = "Based on the configuration of the polymer and state of the last vertex added (occupied or not), and fraction of occupied sites which you get from GrabSites we wil compute the list of weights"

ListOfWeightsCompExp::usage = "Returns a list of weights from the RAWRBFactor_x_x.txt file.  Remember to replace *^ with e in the file BEFORE using this function"

IsOccupied::usage = "Returns whether the site in question is occupied by the polymer.  If it returns {}, site is not occupied by the polymer, if it returns 1 it is and the site is occupied with a protein, if it returns 0 then it is not occupied with a proteine"

OccUnOccPosns::usage = "Returns positions of occupied and unoccupied sites from PolyCoors"

Move::usage = "Returns move based on the index i (1-8) which samples eight possible moves 2*4moves: {{0,-1},{0,1},{-1,0},{1,0},{0,-1},{0,1},{-1,0},{1,0}};"

Rg::usage = "Returns Rg from PolyCoors"

Jump::usage = "Function for protein transfer betweeen two RNAs dressed with protein.  This function will randomly pick an occupied site on one RNA and an unoccupied site on the other, then determine whether the transfer occurs following the Metropolis rule.  Finally it returns the updated polymer coordinates (which includes which sites are occupied with protein)"

Begin["`Private`"]

Jump[polyconfigurationRNA1_,polyconfigurationRNA2_,energy_,epsilon_]:= (

PolyCoorsRNA1=polyconfigurationRNA1;
PolyCoorsRNA2=polyconfigurationRNA2;
Energy=energy;

(*get occupied and unoccupied sites from both RNAs*)
(*start with RNA 1*)
OccUnOccListRNA1=OccUnOccPosns[PolyCoorsRNA1];
(*list of positions for the occupied sites on RNA1*)
OccListRNA1=OccUnOccListRNA1[[1]];
(*number of occupied sites on RNA1*)
numOccRNA1=Length[OccListRNA1];
(*list of positions for the unoccupied sites on RNA1*)
UnOccListRNA1=OccUnOccListRNA1[[2]];

(*now checkout RNA 2*)
OccUnOccListRNA2=OccUnOccPosns[PolyCoorsRNA2];
(*list of positions for the occupied sites on RNA2*)
OccListRNA2=OccUnOccListRNA2[[1]];
(*number of occupied sites on RNA2*)
numOccRNA2=Length[OccListRNA2];
(*list of positions for the unoccupied sites on RNA2*)
UnOccListRNA2=OccUnOccListRNA2[[2]];

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
l=RandomReal[];

focc=N[numOccRNA1/(numOccRNA1+numOccRNA2)];

If[ l <= focc,
      (*attempt to transfer from RNA1 to RNA2*)
  
	(*now we have the positions*)
	If[OccListRNA1 != {} && UnOccListRNA2 != {},
       		 (*so RNA1 has occupied sites and RNA2 has unoccupied sites*)
        	(*pick a random position of an occupied site on RNA1*)
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

			   ,
                           (*reject move*)
			   Energy = Energy;
                          ];

		 ];
  	];
,
	If[OccListRNA2 !={} && UnOccListRNA1 != {},
        	(*so RNA2 has occupied sites and RNA1 has unoccupied sites*)
        	(*pick a random position of an occupied site on RNA2*)
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

			   	,
                           		(*reject move*)
			   		Energy = Energy;
                          	];

		 ];
     ];

];

Return[{PolyCoorsRNA1,PolyCoorsRNA2,Energy}];

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
Return[N[Sqrt[RadGy]/ChainLength]];

)

ListOfWeights[lastvertcoordinate_,polycoordinates_,f_,epsilon_]:= (
focc=f;
PolyCoors=polycoordinates;
LastVertCoor=lastvertcoordinate;
(*let's start sniffing around*)
(*first let's put 0s in the list corresponding to moves that would cause the polymer to run into itself*) 
cycle=0;
FreePosns={};
(*list of weights, constant array of 0s*)
LOW=ConstantArray[0,8];
Do[
TempCoor=LastVertCoor;
TempCoor[[2]]=TempCoor[[2]]+j;
cycle=cycle+1;

OccupancyState=IsOccupied[PolyCoors,TempCoor];
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
TempCoor=LastVertCoor;
TempCoor[[1]]=TempCoor[[1]]+i;
cycle=cycle+1;
(*moving up and down *)
OccupancyState=IsOccupied[PolyCoors,TempCoor];
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
TempCoor=LastVertCoor;

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
	OccupancyState=IsOccupied[PolyCoors,TempCoor2];
	If[OccupancyState == 1,
           numocc=numocc+1;
          ];	

	,{j,{-1,1}}];

	Do[
	TempCoor2=TempCoor;
        TempCoor2[[1]]=TempCoor2[[1]]+i;
	OccupancyState=IsOccupied[PolyCoors,TempCoor2];

        If[OccupancyState == 1,
           numocc=numocc+1;
          ];

         ,{i,{-1,1}}];

Energy=epsilon*numocc;
A=Exp[-Energy];
(*now we need to know how the fraction of occupied sites in the container*)
Weight=A*focc;
(*update list of weights*)
LOW[[Posn]]=Weight;
,

(*these positions correspond to adding an unoccupied site*)

LOW[[Posn]]=(1-focc);

];

,{Posn,FreePosns}];
        
Return[LOW];

)


ListOfWeightsCompExp[rawrbfactorfile_]:= (
(*import data from the RBFactor file*)
LOW=Import[rawrbfactorfile,"List"];
Return[LOW];
)



IsOccupied[polycoordinates_,coordinate_]:= (
PolyCoors=polycoordinates;
Coor=coordinate;
(*get position of the coordinate in the list PolyCoors*)
pos=Position[PolyCoors,Coor];
If[pos == {},
     (*site is not occupied by the polymer*)
     OccupancyState=pos;
     Return[OccupancyState];
  ];

(*change position second element to 2 to get at occupancy state*)
pos[[1,2]]=2;
OccupancyState=PolyCoors[[pos[[1]][[1]],pos[[1]][[2]]]];
Return[OccupancyState];

)

OccUnOccPosns[polymercoordinates_] := (

PolyCoors=polymercoordinates;
OccPosns={};
UnOccPosns={};

numsites=Length[PolyCoors];
Do[
OccupancyState=PolyCoors[[i,2]];
(*Coordinate=PolyCoors[[i,1]];*)

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
Pos=position;
If[Pos < 5,
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
LOW=listofweights;
r=RandomReal[];
(*now implement the select subroutine*)
(*first compute the random cutoff*)
CutOff=Total[N[LOW]]*r;
CumeZ=LOW[[1]];
(*keep increasing CumeZ until we are passed the cutoff*)
i=1;
While[CumeZ < CutOff,
i=i+1;
CumeZ=CumeZ+LOW[[i]];
];
(*return the value i*)
Return[i];
(*This is the position.Make sure you update the distribution of occupied and unoccupied sites at some point, use GrabSites*)

)

Move[position_] := (
Pos=position;
MovesList={{0,-1},{0,1},{-1,0},{1,0},{0,-1},{0,1},{-1,0},{1,0}};

TempMove=MovesList[[Pos]];

Return[TempMove];

)



Grow[polycoordinates_,position_,lastvertcoordinate_,listofweights_,w_,f_]:= (
(*we're assuming you can move, we'll return the updated configuration along with the rosenbluth factor*)
focc=f;
PolyCoors=polycoordinates;
Posn=position;
RBFactor=w;
LastVertCoor=lastvertcoordinate;
LOW=listofweights;
(*list of moves based on position*)


LocalZ=Total[N[LOW]];

(*updated rosenbluth factor, put in the factor of 4 thing*)
RBFactor=(LocalZ*RBFactor)/4;
(*determine the next coordinate*)
 AcceptedMove=Move[Posn];
 NextVertCoor = {LastVertCoor[[1]]+AcceptedMove[[1]],LastVertCoor[[2]]+AcceptedMove[[2]]};

If[Posn < 5,
  
(*site is occupied*)
 NextPolyCoors = {NextVertCoor, 1};
 DelE=-Log[N[LOW[[Posn]]/focc]];
,

 (*site is unoccupied*)
 NextPolyCoors = {NextVertCoor, 0};
 DelE=-Log[N[LOW[[Posn]]/(1-focc)]];
 ];
PolyCoors=Append[PolyCoors,NextPolyCoors];
(*return the updated rosenbluth factor and the updated polymer configuration*)
Return[{PolyCoors, NextVertCoor,RBFactor,DelE}];

)



End[ ]

EndPackage[ ]
