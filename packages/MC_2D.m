BeginPackage["MonteCarlo`"]



MC2DJump::usage = "MC governed jumps between two lattices.  The inputs for this function comes from MC2DTranslate.  One also takes the energy as input."

MC2DAddParticle::usage = "Add a particle to the 2D lattice and return new energy based on nearest-neighbor interactions"

MC2DTranslate::usage = "MC simulation on 2D lattice. Translations only.   Prepare RNA matrix with defective sites.  0-unoccupied, 1-occupied, 2-defective.  As input, give the matrix rnamat, a list of the occupied sites (i.e. {{2,2},{1,4},etc}), the number of steps in the MC, the nearest-neighbor interaction energy epsilon, and the total energy"




Begin["`Private`"]

MC2DJump[rnamat1_,rnamat2_,epsilon_,energy_] := (

(*make local variables*)
RNA1mat=rnamat1;
RNA2mat=rnamat2;
Energy=energy;

NumRowsRNA1=Dimensions[RNA1mat][[1]];
NumColumnsRNA1=Dimensions[RNA1mat][[2]];

NumRowsRNA2=Dimensions[RNA2mat][[1]];
NumColumnsRNA2=Dimensions[RNA2mat][[2]];

OccupiedSitesRNA1=Position[RNA1mat,1];
UnOccupiedSitesRNA2=Position[RNA2mat,0];

(*pick random occupied site on RNA1 and unoccupiedsite RNA2*)
RandomOccRNA1=OccupiedSitesRNA1[[RandomInteger[{1,Length[OccupiedSitesRNA1]}]]];
RandomUnOccRNA2=UnOccupiedSitesRNA2[[RandomInteger[{1,Length[UnOccupiedSitesRNA2]}]]];



If[ RandomOccRNA1 == {} || RandomUnOccRNA2 == {},
       (*can't move*)
       Return[{RNA1mat,RNA2mat,Energy}];
  ];

(*Calculate initial energy*)
(*we must determine the occupied and unoccupied nearest-neighbors*)
NNOccSitesRNA1=0;


(*RNA*)

	Do[
		If[RNA1mat[[Mod[RandomOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomOccRNA1[[2]]+j,NumColumnsRNA1,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomOccRNA1[[2]]+j,NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA1mat[[Mod[RandomOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomOccRNA1[[2]],NumColumnsRNA1,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomOccRNA1[[2]],NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

InitialEnergy=NNOccSitesRNA1*epsilon;




NNOccSitesRNA2=0;
	Do[
		If[RNA2mat[[Mod[RandomUnOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]]+j,NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomUnOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]]+j,NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA2mat[[Mod[RandomUnOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]],NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomUnOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomUnOccRNA2[[2]],NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

FinalEnergy=NNOccSitesRNA2*epsilon;

DelE=FinalEnergy-InitialEnergy;

(*MC time!*)
If[DelE <= 0,
         (*accept move*)
	 RNA2mat[[RandomUnOccRNA2[[1]],RandomUnOccRNA2[[2]]]]=1;
         RNA1mat[[RandomOccRNA1[[1]],RandomOccRNA1[[2]]]]=0;
	Energy=Energy+DelE;
	,
	(*accept move with boltzman probability*)
       	r = RandomReal[];
	If[r<Exp[-(DelE)],
             (*accept move*)	
	 	RNA2mat[[RandomUnOccRNA2[[1]],RandomUnOccRNA2[[2]]]]=1;
        	RNA1mat[[RandomOccRNA1[[1]],RandomOccRNA1[[2]]]]=0;
		Energy=Energy+DelE;
		,
		(*reject move*)
		Energy=Energy;
	];

];


OccupiedSitesRNA2=Position[RNA2mat,1];
UnOccupiedSitesRNA1=Position[RNA1mat,0];

(*pick random occupied site on RNA2 and unoccupiedsite RNA1*)
RandomOccRNA2=OccupiedSitesRNA2[[RandomInteger[{1,Length[OccupiedSitesRNA2]}]]];
RandomUnOccRNA1=UnOccupiedSitesRNA1[[RandomInteger[{1,Length[UnOccupiedSitesRNA1]}]]];

(*pick random occupied site on RNA2 and unoccupied site on RNA1*)
If[ RandomOccRNA2 == {} || RandomUnOccRNA1 == {},
       (*can't move*)
       Return[{RNA1mat,RNA2mat,Energy}];
  ];

(*Calculate initial energy*)
(*we must determine the occupied and unoccupied nearest-neighbors*)
NNOccSitesRNA2=0;


(*RNA*)

	Do[
		If[RNA2mat[[Mod[RandomOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomOccRNA2[[2]]+j,NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomOccRNA2[[1]],NumRowsRNA2,1],Mod[RandomOccRNA2[[2]]+j,NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA2mat[[Mod[RandomOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomOccRNA2[[2]],NumColumnsRNA2,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA2mat[[Mod[RandomOccRNA2[[1]]+i,NumRowsRNA2,1],Mod[RandomOccRNA2[[2]],NumColumnsRNA2,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA2 = NNOccSitesRNA2+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

InitialEnergy=NNOccSitesRNA2*epsilon;



NNOccSitesRNA1=0;
	Do[
		If[RNA1mat[[Mod[RandomUnOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]]+j,NumColumnsRNA1,1]]] != 2,
		(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomUnOccRNA1[[1]],NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]]+j,NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNA1mat[[Mod[RandomUnOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]],NumColumnsRNA1,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNA1mat[[Mod[RandomUnOccRNA1[[1]]+i,NumRowsRNA1,1],Mod[RandomUnOccRNA1[[2]],NumColumnsRNA1,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSitesRNA1 = NNOccSitesRNA1+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

FinalEnergy=NNOccSitesRNA1*epsilon;

DelE=FinalEnergy-InitialEnergy;

(*MC time!*)
If[DelE <= 0,
         (*accept move*)
	 RNA2mat[[RandomOccRNA2[[1]],RandomOccRNA2[[2]]]]=0;
         RNA1mat[[RandomUnOccRNA1[[1]],RandomUnOccRNA1[[2]]]]=1;
	Energy=Energy+DelE;
	,
	(*accept move with boltzman probability*)
       	r = RandomReal[];
	If[r<Exp[-(DelE)],
             (*accept move*)	
	 	RNA2mat[[RandomOccRNA2[[1]],RandomOccRNA2[[2]]]]=0;
        	RNA1mat[[RandomUnOccRNA1[[1]],RandomUnOccRNA1[[2]]]]=1;
		Energy=Energy+DelE;
		,
		(*reject move*)
		Energy=Energy;
	];

];
 

Join[{RNA1mat},{RNA2mat},{Energy}]






)








MC2DAddParticle[rnamat_,epsilon_,energy_] := (
(*make local copies of variables*)
RNAmat=rnamat;
Energy=energy;
NumRows=Dimensions[RNAmat][[1]];
NumColumns=Dimensions[RNAmat][[2]];

(*Find unoccupiedsites*)

UnOccupiedSites=Position[RNAmat,0];

If[UnOccupiedSites == {},
   Return[{RNAmat,Energy}];
  ];

RandomUnOcc=UnOccupiedSites[[RandomInteger[{1,Length[UnOccupiedSites]}]]];

DelE=0;
(*now calculate the number of occupied nearest-neighbors*)
NNOccSites=0;

	Do[
		If[RNAmat[[Mod[RandomUnOcc[[1]],NumRows,1],Mod[RandomUnOcc[[2]]+j,NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomUnOcc[[1]],NumRows,1],Mod[RandomUnOcc[[2]]+j,NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = NNOccSites+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNAmat[[Mod[RandomUnOcc[[1]]+i,NumRows,1],Mod[RandomUnOcc[[2]],NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomUnOcc[[1]]+i,NumRows,1],Mod[RandomUnOcc[[2]],NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = NNOccSites+1;
					,
					(*nearest-neighbor site is unoccupied*)
		          		Energy=Energy;	
			];

		];

	,{i,-1,1,2}];

DelE=NNOccSites*epsilon;
Energy=Energy+DelE;
(*update matrix*)
RNAmat[[RandomUnOcc[[1]],RandomUnOcc[[2]]]]=1;

Join[{RNAmat},{Energy}]

)












MC2DTranslate[rnamat_, numsteps_,epsilon_,energy_] := (


(*make local variables of RNAmat,occupiedsitesRNA,occupiedsitesRNA2*)
RNAmat=rnamat;
NumRows=Dimensions[RNAmat][[1]];
NumColumns=Dimensions[RNAmat][[2]];
Energy=energy;

OccupiedSites=Position[RNAmat,1];

UnOccupiedSites=Position[RNAmat,0];

If[OccupiedSites == {} || UnOccupiedSites == {},
   Return[{RNAmat,Energy}];];

(*LOOP TIME BABY*)

Do[
(*pick random occupied sites on RNA *)

RandomOcc=OccupiedSites[[RandomInteger[{1,Length[OccupiedSites]}]]];

(*we must determine the occupied and unoccupied nearest-neighbors*)
NNOccSites={};
NNUnOccSites={};

(*RNA*)

	Do[
		If[RNAmat[[Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = Append[NNOccSites,{Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		NNUnOccSites = Append[NNUnOccSites,{Mod[RandomOcc[[1]],NumRows,1],Mod[RandomOcc[[2]]+j,NumColumns,1]}];	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNAmat[[Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					NNOccSites = Append[NNOccSites,{Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		NNUnOccSites = Append[NNUnOccSites,{Mod[RandomOcc[[1]]+i,NumRows,1],Mod[RandomOcc[[2]],NumColumns,1]}];	
			];

		];

	,{i,-1,1,2}];
(*initial energy*)
InitialEnergy=Length[NNOccSites]*epsilon;

(*pick random unoccupied nearest-neighbor*)
RandomNNUnOccSites=NNUnOccSites[[RandomInteger[{1,Length[NNUnOccSites]}]]];

If[RandomNNUnOccSites == {},
       (*can't move*)
	FinalEnergy=InitialEnergy;
	OccNNNN={};
			,
	
	(*not empty, so make a virtual move*)	
	
         (*figure out the occupied nearest-neighbors of the randomly chosen unoccupied nearest neighbor NOT including RandOcc*)
	OccNNNN={};
	
	Do[
		If[RNAmat[[Mod[RandomNNUnOccSites[[1]],NumRows,1],Mod[RandomNNUnOccSites[[2]]+j,NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomNNUnOccSites[[1]],NumRows,1],Mod[RandomNNUnOccSites[[2]]+j,NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					OccNNNN = Append[OccNNNN,{Mod[RandomNNUnOccSites[[1]],NumRows,1],Mod[RandomNNUnOccSites[[2]]+j,NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		(*we don't need to keep track of them*)
					Energy=Energy;	
			];

		];

	,{j,-1,1,2}];

	Do[
		If[RNAmat[[Mod[RandomNNUnOccSites[[1]]+i,NumRows,1],Mod[RandomNNUnOccSites[[2]],NumColumns,1]]] != 2,
			(*nearest-neighbor site is NOT defective*)
			If[RNAmat[[Mod[RandomNNUnOccSites[[1]]+i,NumRows,1],Mod[RandomNNUnOccSites[[2]],NumColumns,1]]] == 1,
					(*nearest-neighbor site is occupied*)
					OccNNNN = Append[OccNNNN,{Mod[RandomNNUnOccSites[[1]]+i,NumRows,1],Mod[RandomNNUnOccSites[[2]],NumColumns,1]}];
					,
					(*nearest-neighbor site is unoccupied*)
		          		(*we don't need to keep track of them*)	
					Energy=Energy;
			];

		];

	,{i,-1,1,2}];



	(*but this procedure will have included RandomOcc*)

        
	FinalEnergy=Length[Complement[OccNNNN,{RandomOcc}]]*epsilon;

  ];

(*Now calculate change in Energy*)
DelE=FinalEnergy-InitialEnergy;

(*Now make monte-carlo move*)
If[DelE <= 0,
   (*accept move*)
  	  If[RandomNNUnOccSites != {},
	    (*we can move*)
	    Energy=Energy+DelE;
	    (*swap occupied site for the randomly chosen unoccupied site*)
	     OccupiedSites=Complement[Append[OccupiedSites,RandomNNUnOccSites],{RandomOcc}];
	     UnOccupiedSites=Complement[Append[UnOccupiedSites,RandomOcc],{RandomNNUnOccSites}];
	(*update matrix*)
	    RNAmat[[RandomNNUnOccSites[[1]],RandomNNUnOccSites[[2]]]]=1;
	    RNAmat[[RandomOcc[[1]],RandomOcc[[2]]]]=0;
		,
	     (*can't move*)
	     Energy=Energy;
	     OccupiedSites=OccupiedSites;
	     UnOccupiedSites=UnOccupiedSites;
	    ];
	,

	(*DelE is > 0 so accept with boltzman proabability*)
	r = RandomReal[];
	If[r<Exp[-(DelE)],
		(*accept move*)
		Energy=Energy+DelE;
		OccupiedSites=Complement[Append[OccupiedSites,RandomNNUnOccSites],{RandomOcc}];
	        UnOccupiedSites=Complement[Append[UnOccupiedSites,RandomOcc],{RandomNNUnOccSites}];
	    (*update matrix*)
	    RNAmat[[RandomNNUnOccSites[[1]],RandomNNUnOccSites[[2]]]]=1;
	    RNAmat[[RandomOcc[[1]],RandomOcc[[2]]]]=0;
		,
		(*reject move*)
		Energy=Energy;
		OccupiedSites=OccupiedSites;
		UnOccupiedSites=UnOccupiedSites;


          ];

 ];

,{i,1,numsteps}];

Join[{RNAmat},{Energy}]

)

End[ ]

EndPackage[ ]
