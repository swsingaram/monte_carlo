(* ::Package:: *)

Off[General::compat]
<< Combinatorica`; << GraphUtilities`;
RTree[n_] := (System`AdjacencyGraph[
    Combinatorica`ToAdjacencyMatrix[Combinatorica`RandomTree[n]]]);
PShuffle[t_] := (System`AdjacencyGraph[
    Combinatorica`ToAdjacencyMatrix[
     Combinatorica`CodeToLabeledTree[
      System`RandomSample[
       Combinatorica`LabeledTreeToCode[
        GraphUtilities`ToCombinatoricaGraph[
         System`AdjacencyMatrix[t]]]]]]]);
Dendrimer[f_, g_] := 
 Module[{tree, vtree, subtree, jsite, fulltree}, 
  tree = System`CompleteKaryTree[g, f - 1]; 
  vtree = System`VertexCount[tree]; jsite = vtree + 1; 
  subtree = System`CompleteKaryTree[g - 1, f - 1]; 
  fulltree = System`GraphDisjointUnion[tree, subtree];
  fulltree = System`EdgeAdd[fulltree, UndirectedEdge[1, jsite]]]
Rg[t_] := 
  Module[{eli, eln, prodl, tl, ed, new, ll, rr, pr, tprod}, 
   eli = System`EdgeList[t]; eln = Length[eli]; prodl = {};
   tl = VertexCount[t];
   Do[ed = Take[eli, {j}]; new = EdgeDelete[t, ed]; 
    {ll, rr} = System`ConnectedComponents[new];
    pr = Length[ll]*Length[rr]; 
    prodl = Append[prodl, pr]; 
    , {j, eln}
    ];
   tprod = Total[prodl]; 
   N[Sqrt[tprod/(tl^2)]]
   ];
PPairShuffle[t_, s_] := 
 Module[{bp, nmax, picks, pick1, pick2}, 
  bp = Combinatorica`LabeledTreeToCode[
    GraphUtilities`ToCombinatoricaGraph[System`AdjacencyMatrix[t]]];
  nmax = Length[bp];
  Do[
   picks = RandomChoice[Range[nmax], 2];
   pick1 = bp[[picks[[1]]]];
   pick2 = bp[[picks[[2]]]];
   bp[[picks[[2]]]] = pick1;
   bp[[picks[[1]]]] = pick2;
   , {s}];
  System`AdjacencyGraph[
   Combinatorica`ToAdjacencyMatrix[
    Combinatorica`CodeToLabeledTree[bp]]]
  ]
RgEnsemble[t_, n_] := Module[{dat, r, em, es, z},
  dat = {}; Do[dat = Append[dat, Rg[PShuffle[t]]];, {n}];
  r = Rg[t]; em = Mean[dat]; es = StandardDeviation[dat]; 
  z = (r - em)/es;
  {VertexCount[t], r, em, es, z}
  ]
ErodeTree[t_] := 
 Module[{a, aa, av, ad, adp, n, ag, tab, dat, ae, vc, cycle, a1p},
  n = VertexCount[t];
  dat = {}; tab = {};
  ae = t;
  cycle = 0;
  tab = Append[tab, {0, n}];
  vc = n;
  While[vc > 2,
   cycle = cycle + 1;
   av = System`VertexList[ae];
   ad = VertexDegree[ae];
   {adp} = Transpose[Position[ad, 1]];
   a1p = av[[adp]];
   ae = VertexDelete[ae, a1p];
   vc = VertexCount[ae];
   dat = {cycle, vc};
   tab = Append[tab, dat];
   
   ];
  tab
  ]
ErodeLeaves[t_] := 
 Module[{a, aa, av, ad, adp, n, ag, tab, dat, ae, vc, cycle, a1p},
  n = VertexCount[t];
  dat = {}; tab = {};
  ae = t;
  cycle = 0;
  vc = n;
  While[vc > 2,
   cycle = cycle + 1;
   av = System`VertexList[ae];
   ad = VertexDegree[ae];
   {adp} = Transpose[Position[ad, 1]];
   a1p = av[[adp]];
   ae = VertexDelete[ae, a1p];
   dat = {cycle, Length[a1p]};
   tab = Append[tab, dat];
   vc = VertexCount[ae];
   ];
  tab
  ]

(* REQUIRED SUB MODULE*)

LoopPairExtract[refpair_, hstarts_] := 
  Module[{looppairs, refstart, refend, counter, pairpos, pair, begin, 
    post, hs},
   looppairs = {};
   hs = hstarts;
   looppairs = Append[looppairs, refpair];
   refstart = refpair[[1]];
   refend = refpair[[2]];
   counter = refstart + 1;
   For[counter == refstart + 1, counter < refend, counter++,
    If[Length[Position[hs, counter]] == 1,
      {{pairpos, post}} = Position[hs, counter];
      pair = hs[[pairpos]];
      hs = Delete[hs, pairpos];
      looppairs = Append[looppairs, pair];
      counter = pair[[2]];
      ];
    ];
   {looppairs, hs}];

(* MAIN MODULE*)

CT2Tree[ff_] := 
 Module[{cta, hstarts, hends, ctb, tree, lp, loops, seql, hstartscopy,
    hendscopy, refpair, tolist, vto, e, f, dlist, ct, b},
  
  (*OPEN CT FILE AND READ BASE-PAIRINGS*)
  If[FileExistsQ[ff],
   ct = Import[ff, "Table"], Print["File " <> ff <> " does not exist"];
   ];
  cta = ct[[All, 1]];
  ctb = ct[[All, 5]];
  
  (* REMOVE CT FILE HEADER LINE *)
  cta = Delete[cta, 1];
  ctb = Delete[ctb, 1];
  
  (* PICK UNIQUE BASEPAIRS AT START AND END OF HELICES *)
  hends = {};
  hstarts = {};
  seql = Length[cta];
  Do[
   If[i == 1  && ctb[[i]] != 0 &&  ctb[[i + 1]] == ctb[[i]] - 1, 
    hstarts = Append[hstarts, {cta[[i]], ctb[[i]]}]
    ];
   If[ctb[[i]] != 0 && 
     1 < i < seql && (ctb[[i - 1]] == 0 || 
       ctb[[i - 1]] != ctb[[i]] + 1) && ctb[[i + 1]] == ctb[[i]] - 1  && 
      cta[[i]] <=  ctb[[i]],
    hstarts = Append[hstarts, {cta[[i]], ctb[[i]]}] 
    ];
   If[ctb[[i]] != 0 && 
     1 < i < seql && (ctb[[i + 1]] == 0 || 
       ctb[[i + 1]] != ctb[[i]] - 1 ) && 
     ctb[[i - 1]] == ctb[[i]] + 1 &&  cta[[i]] <=  ctb[[i]],
    hends = Append[hends, {cta[[i]], ctb[[i]]}] 
    ];
   , {i, 1, seql}
   ];(*end do over i*)
  
  (* REMOVE HELICES SHORTER THAN 3  
  dlist={};
  Do[
	e=hends[[l,1]];
	f=hstarts[[l,1]];
	If[e-f < 2,
	dlist=Append[dlist,{l}];
	];
	,{l,1,Length[hends]}
	];
  hstarts=Delete[hstarts,dlist];
  hends=Delete[hends,dlist];*)
  hstartscopy = hstarts;
  hendscopy = hends;

  (* CHECK FOR EQUAL NUMBER OF STARTS AND ENDS *)
  
  If[Length[hstarts] != Length[hends], 
   Print["Unequal Half-Edges - Check CT file. Exiting."]; Exit
   ];
  
  
  (* MAKE LIST OF INCOMING AND OUTGOING HELICES AT LOOPS *)
  
  loops = {};
  
  (*External Loop*)
  refpair = {0, seql};
  lp = Part[LoopPairExtract[refpair, hstartscopy], 1];
  (*lp=Delete[lp,1];*)
  loops = Append[loops, lp];
  
  (*Internal Loops*)
  While[Length[hendscopy] >= 1,
   refpair = hendscopy[[1]];
   {lp, hstartscopy} = LoopPairExtract[refpair, hstartscopy];
   hendscopy = Delete[hendscopy, 1];
   loops = Append[loops, lp];
   ];
  
  (* Remove loops with less than 2 bp (consistent with RAG)
 LoopLength[lp_]:= Module[{length,a,b},
 Do[
    length=0;
   If[ll==1, 
	a=lp[[{1,1}]]; b=lp[[{2,1}]];
    length=length+(b-a);
	];
	If[ll>1,
	a=lp[[{ll,2}]]; b=lp[[{ll+1,1}]];
    length=length+(b-a);
	];
  ,{ll,Length[lp]}];
  length
  ];
  
  loopsizes=Map[LoopLength[#]&,loops];
  sbb=Position[loopsizes,1];*)

  (* CREATE TREE *)
  tree = {};
  b = {};
  Do[
   tolist = Map[Position[hstarts, #] + 1 &, Delete[loops[[v]], 1]];
   Do[{{vto}} = tolist[[t]];
    tree = Append[tree, UndirectedEdge[v, vto]]
    , {t, Length[tolist]}
    ];
   , {v, 1, Length[hends]}
   ];
  System`Graph[tree, GraphLayout -> "SpringEmbedding"]
  ];

ProxPair[ff_] := 
   Module[{ct, cta, pprox, ctb, seql},
   (*OPEN CT FILE AND READ BASE-PAIRINGS*)
   If[FileExistsQ[ff],
    ct = Import[ff, "Table"], 
    Print["File " <> ff <> " does not exist"];
    ];
   cta = ct[[All, 1]];
   ctb = ct[[All, 5]];
   (* REMOVE CT FILE HEADER LINE *)
   cta = Delete[cta, 1];
   ctb = Delete[ctb, 1];
   pprox = {};
   seql = Length[cta];
   Do[
    If[cta[[i]] < ctb[[i]],
      pprox = Append[pprox, (ctb[[i]] - cta[[i]])]
      ];
    , {i, 1, seql}
    ];
   {N[Mean[pprox]], N[StandardDeviation[pprox]]}
   ];
