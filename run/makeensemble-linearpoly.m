Needs["MonteCarloRosenbluthBranchedPolymer`","~/scripts/mathematica/packages/MC_Rosenbluth-branchedpoly.m"];
(*let us  grow a linear polymer without any bound sites*)
$HistoryLength=1;
(*number of attempted generations*)
i=1;
numruns=1;
While[i <= numruns,

(*specify the number of NumOcc and NumUnOcc initially and energy*)
NumOcc=10;
NumUnOcc=10;
epsilon=-3;

f=N[NumOcc/(NumOcc+NumUnOcc)];
CPRNAratio=f;
(*generate a number between 0 and 1*)
r=RandomReal[];
Energy=0;

If[r < f,
    (*randomly selected an occupied site*)
    grablist=GrabSites[NumOcc,NumUnOcc,1];
    NumOcc=grablist[[1]];
    NumUnOcc=grablist[[2]];
    (*initialize coordinates*)
    PolyCoors={{{0,0},1,1}};
   f=N[NumOcc/(NumOcc+NumUnOcc)];
     ,
    (*randomly selected an unoccupied site*)
    grablist=GrabSites[NumOcc,NumUnOcc,5];
    NumOcc=grablist[[1]];
    NumUnOcc=grablist[[2]];
    (*initialize coordinates*)
    PolyCoors={{{0,0},0,1}};
   f=N[NumOcc/(NumOcc+NumUnOcc)];
  ];

Steps=NumOcc+NumUnOcc;
LastVertCoor={0,0};
RBFactor=1;
vertindex=1;
Do[
vertindex=vertindex+1;
(*calculate list of weights*)
LOW=ListOfWeights[LastVertCoor,PolyCoors,f,epsilon];
(*check if LOW is all 0*)
 If[Total[N[LOW]] == 0,
          (*got stuck*)
          (*append to file empty polycoors*)
          PolyCoors={};
          RBFactor=0;
          Return[{PolyCoors,RBFactor}];
   ];
(*now select a configuration*)
PosnList=ConfigSelect[LOW];
(*now we need to update NumOcc,NumUnOcc, and Grow*)
GrowList=Grow[PolyCoors,PosnList,vertindex,LastVertCoor,LOW,RBFactor,f];
(*update polymer*)
PolyCoors=GrowList[[1]];
(*next or frontier coordinates*)
LastVertCoor=GrowList[[2]];
RBFactor=GrowList[[3]];
(*update energy*)
DelE=GrowList[[4]];
Energy=Energy+DelE;
(*update NumOcc,NumUnOcc*)
grablist=GrabSites[NumOcc,NumUnOcc,Posn];
NumOcc=grablist[[1]];
NumUnOcc=grablist[[2]];
(*update fraction of occupied sites*)
If [ (NumOcc+NumUnOcc) != 0,
f=N[NumOcc/(NumOcc+NumUnOcc)];
   ];

(*clear variables*)
ClearAll[grablist,Growlist,grablist,LOW,PosnList];

,{stepnum,1,Steps}];
(*return to write PolyCoors to file*)
If[RBFactor != 0,
PutAppend[PolyCoors,StringJoin["RAWPolyCoors_",ToString[CPRNAratio],"_",ToString[Steps],".txt"]];
PutAppend[Rg[PolyCoors],StringJoin["radiusgyration_",ToString[CPRNAratio],"_",ToString[Steps],".txt"]];
PutAppend[RBFactor,StringJoin["RAWRBFactor_",ToString[CPRNAratio],"_",ToString[Steps],".txt"]];
PutAppend[Energy,StringJoin["RAWEnergy_",ToString[CPRNAratio],"_",ToString[Steps],".txt"]];

(*increment i*)
i=i+1;

  ];
(*clear variables*)
ClearSystemCache[];

(*end of do loop*)
];
Exit[]
