<<~/scripts/mathematica/packages/protdist.m
Do[
   Sigma09280 = ProtDistCalc[0.9,300,i,1000,280];
   Sigma06280 = ProtDistCalc[0.6,300,i,1000,280];
   Sigma06180 = ProtDistCalc[0.6,300,i,1000,180];
   
   Sigma09280 >>> 09280_protdist.dat;
   Sigma06280 >>> 06280_protdist.dat;
   Sigma06180 >>> 06180_protdist.dat;
   , {i,1,1000}];

Do[
   ProtTit = ProtDistCalc[j,300,800,1000,180];
   ProtTit >>> ProteinTitration_80c_180.dat;
  , {j,0.1,0.9,0.1}];

Exit[]
