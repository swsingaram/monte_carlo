BeginPackage["MonteCarlo`"]




MCEnsembleOneD::usage = "MCEnsembleOneD[m,n,no,trials,step,e] runs the Metropolis Algorithm on m 1D lattices each with n sites.  The total # of proteins in the system is no.  trials is the number of trials and step is the number of trials between two samples, e is the lateral energy"

 Begin["`Private`"]

MCEnsembleOneD[m_,n_,no_,trials_,step_,e_] := (OccupiedSites = {}; While[Less[Length[OccupiedSites], no],OccupiedSites = Union[OccupiedSites,DeleteDuplicates[RandomSample[Range[Times[n,m]], no]]];UnOccupiedSites = Complement[Range[Times[n,m]], OccupiedSites];Config = ConstantArray[0,Times[n,m]];Energy = 0;Do[Config[[OccupiedSites[[i]]]] = 1;If[And[Unequal[Mod[OccupiedSites[[i]],n],1], Equal[Config[[OccupiedSites[[i]] - 1]], 1]],Energy = Energy - e;];,{i,1,no}];Do[TempConfig = Config;TempOcc = OccupiedSites;  TempUnOcc = UnOccupiedSites; randI1 = RandomInteger[{1,no}];  randOcc = OccupiedSites[[randI1]];  randI2 = RandomInteger[{1, Times[m,n] - no}];  randUnOcc = UnOccupiedSites[[randI2]];  TempConfig[[randOcc]] = 0;  TempConfig[[randUnOcc]] = 1; TempOcc[[randI1]] = randUnOcc;  TempUnOcc[[randI2]] = randOcc;  Dele = 0;  nni = 0;  nnf = 0;    If[And[Unequal[Mod[randOcc,n],0],Equal[Config[[randOcc + 1]], 1]],  nni = nni + 1;,  nni = nni;  ];    If[And[Unequal[Mod[randOcc,n],1],Equal[Config[[randOcc - 1]], 1]],      nni = nni + 1;,      nni = nni;      ];    If[And[Unequal[Mod[randUnOcc,n],0],Equal[TempConfig[[randUnOcc + 1]], 1]],      nnf = nnf + 1;,      nnf = nnf;     ];  If[And[Unequal[Mod[randUnOcc,n],1],Equal[TempConfig[[randUnOcc - 1]], 1]],  nnf = nnf + 1;,  nnf = nnf;  ];Delnn = nnf - nni;Dele = Times[Delnn,Times[-1,e]];  If[LessEqual[Dele,0],    Energy = Energy + Dele;   Config = TempConfig;    OccupiedSites = TempOcc;    UnOccupiedSites = TempUnOcc;,    r = RandomReal[];    If[Less[r,Exp[Times[-1,Dele]]],       Energy = Energy + Dele;       Config = TempConfig;       OccupiedSites = TempOcc;       UnOccupiedSites = TempUnOcc;,       Energy = Energy;       Config = Config;       OccupiedSites = OccupiedSites;       UnOccupiedSites = UnOccupiedSites;      ];     ];If[Equal[Mod[numtrials,step],0], Do[    WriteConfig = Part[Config,((i - 1)*n + 1);;(i*n)]; NumProt = Count[WriteConfig,1]; PutAppend[WriteConfig,StringJoin[ToString[Divide[no,m]],"configurations.dat"]]; PutAppend[NumProt,StringJoin[ToString[Divide[no,m]],"numprotein.dat"]]; ,{i,1,m}];  PutAppend[Energy,StringJoin[ToString[Divide[no,m]],"energy.dat"]]; ];,{numtrials,1,trials}];];)


End[ ]

EndPackage[ ]
