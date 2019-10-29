Needs["MonteCarlo`","~/scripts/mathematica/packages/MC.m"];
(*ParallelEvaluate[Needs["MonteCarlo`","~/scripts/mathematica/packages/MC.m"];];*)
(*LaunchKernels[4];*)
(*50% sites defective*)
(*epsilonRNA2=-5;*)
(*epsilonRNA1=-3.75;*)
(*numsteps=100;*)
numjumps=100000000;
$HistoryLength=1;
(*DistributeDefinitions[numjumps];*)

rnamat1={{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
Do[
(**)
rnamat2={{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};


Print[Dimensions[rnamat2]];
(*add particles and translate a bit*)
(*we will make  them 50% covered *)

Energy = 0;
Do[
AddedParticleList = MC2DAddParticle[rnamat2,epsilonRNA2,Energy];
(*a particle has been added to rnamat2*)
rnamat2=AddedParticleList[[1]];
Energy=AddedParticleList[[2]];

(*TranslateStepList = MC2DTranslate[rnamat2,numsteps,epsilonRNA2,Energy];
rnamat2=TranslateStepList[[1]];
Energy=TranslateStepList[[2]];*)

,{particlenum,1,75}];

Print["Number of occupied sites on rna 2"];
Print[Length[Position[rnamat2,1]]];
Print["*****************"];
(*ClearAll[MC2DAddParticle];*)

(*

Do[
AddedParticleList = MC2DAddParticle[rnamat1,epsilonRNA1,Energy];
(*a particle has been added to rnamat1*)
rnamat1=AddedParticleList[[1]];
Energy=AddedParticleList[[2]];

TranslateStepList = MC2DTranslate[rnamat1,numsteps,epsilonRNA1,Energy];

rnamat1=TranslateStepList[[1]];
Energy=TranslateStepList[[2]];


,{particlenum,1,37}];

Print["Number of occupied sites on rna 1"];
Print[Length[Position[rnamat1,1]]];
Print["*****************"];

*)

Do[

(*If[ Position[rnamat2,0] == {} || Position[rnamat1,0] == {} || Position[rnamat2,1] == {} || Position[rnamat1,1] == {},
        Break[];]; *)

(*particles added now allow them to translate and jump*)

(*Translate rna2*)
(*TranslateStepList = MC2DTranslate[rnamat2,numsteps,epsilonRNA2,Energy];
rnamat2=TranslateStepList[[1]];
Energy=TranslateStepList[[2]];*)

(*(*Translate rna1*)
TranslateStepList = MC2DTranslate[rnamat1,numsteps,epsilonRNA1,Energy];
rnamat1=TranslateStepList[[1]];
Energy=TranslateStepList[[2]];
*)

(*Now jump*)
JumpStepList = MC2DJump[rnamat2,rnamat2,epsilonRNA2,epsilonRNA2,Energy];
(*rnamat1=JumpStepList[[1]];*)
rnamat2=JumpStepList[[2]];
Energy=JumpStepList[[3]];


(*PutAppend[rnamat1,"rnamat1.dat"];*)
If[ Mod[i,100000] == 0,
PutAppend[rnamat2,StringJoin["rnamat2_",ToString[epsilonRNA2],".dat"]];
PutAppend[Energy,StringJoin["energy_",ToString[epsilonRNA2],".dat"]];
  ];

,{i,1,numjumps}];

Print["ALL DONE!!!"];
(*clear variables and cache*)
ClearSystemCache[];
Clear[Energy,TranslateStepList,JumpStepList,AddedParticleList];
(*clear functions and reload at the beginning of the do loop*)
(*ClearAll[MC2DJump,MC2DTranslate];*)

,{epsilonRNA2,{0,-1,-2,-3,-4}}];

Exit[] 
