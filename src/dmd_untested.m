(* ::Package:: *)

(* ::Title::Closed:: *)
(*untested*)


(* ::Subsubsection::Closed:: *)
(*1/z without breaking the system*)


carefulinverse[list_]:=(list + $MachineEpsilon*Exp[I*Arg[list]])^(-1);


(* ::Subsubsection::Closed:: *)
(*Snapshots from an IVP*)


getimeseries[ic_,vfield_,tinit_,tfinn_,sampackage_,chosenputty_]:=samplethistraj[(x/.Flatten[NDSolve[{x'[t]==vfield[t,x[t]],x[tinit]==ic},x,{t,tinit,tfinn}]]),sampackage,chosenputty];


(* ::Subsubsection::Closed:: *)
(*Combining mode-norms and Jacobians into an eigen-value score*)


wtsnjacs2rfactors[wts_,jacs_,{npow_,dpow_}]:=wts^npow/(1+jacs^dpow);


(* ::Subsection::Closed:: *)
(*Dictionaries*)


(* ::Subsubsection::Closed:: *)
(*Polynomial *)


naiveboxevals[x_,{widths_,numvec_,nboxes_,boxvol_,boxranges_,trueboxranges_}]:=Module[{index},
index = (Floor@((x - boxranges[[All,1]])/widths+1/2)).numvec + 1;
SparseArray[{index-> 1/Sqrt[boxvol]},{nboxes}]
];


(* ::Text:: *)
(*splpower: To take care of the rare possibility that you have 0^0 during pascalift*)


splpower[x_,y_]:=If[y !=0,Power[x,y],0*x+1];
SetAttributes[splpower,Listable];


(* ::Subsubsection::Closed:: *)
(*triglift scaling it before lifting*)


triglift[mat_,{{scalelist_},trueboxranges_}]:=Module[{scaledmat},
scaledmat =((trueboxranges[[All,2]])^-1)*mat;
(* Scale - SinCos - Outer to take all possible combinations - Flatten to make each individual row stand out again *)
Flatten[
Map[
(Flatten[ 
((Outer[Times,##,1]&)@@Map[{Sin[\[Pi]*#],Cos[\[Pi]*#]}&,#]),Length[#]-1]&)[splscale[#,scaledmat]]&,scalelist],1]
];


naivebasislift[mat_,params_,flavour_]:=Which[
(* Delay / Monomials + Pascal's triangle*)
needpolytri[flavour],pascalift[mat,params[[1]]],
(* trigtri business *)
TrueQ[flavour=="trigtri"],triglift[mat,params],
(* All other dictionaries *)
True,Map[lift[#,params[[1]],flavour]&,mat\[Transpose]]\[Transpose]
];


(* ::Subsubsection::Closed:: *)
(*Using insideQnD to judge eval points*)


basislift[mat_,params_,flavour_]:=Module[{logiclist,insideones,liftedones,outlifted,index},
If[flavour=="boxfuns",
(* In/out check is implicit with box functions *)
naivebasislift[mat,params,flavour],
(* Now, need an explicit check *)
logiclist = Map[insideQnD[#,params[[2]]]&,mat\[Transpose]];
insideones = Pick[mat\[Transpose],logiclist]\[Transpose];
liftedones = naivebasislift[insideones,params,flavour]\[Transpose];
outlifted = liftedones[[1]]*0;
index = 1;
((Reap[Scan[If[#,Sow[liftedones[[index]]];index++,Sow[outlifted]]&,logiclist]])[[2,1]])\[Transpose]
]
];


cookboxgrub[{boxranges_,centercounts_}]:=Module[{widths,numvec,nboxes,boxvol,trueboxranges},
widths = Map[(#.{-1,1})&,boxranges]/(centercounts-1);
(* Indexing so taht you count the number of boxes as approrpiate for each coordinate *)
numvec = FoldList[Times,1,Most@centercounts];
(* Self evident parameters *)
nboxes= Times@@centercounts;
boxvol = Times@@widths;
trueboxranges =boxranges + {widths}\[Transpose].{{-1/2,1/2}}; 
{widths,numvec,nboxes,boxvol,boxranges,trueboxranges}
];
