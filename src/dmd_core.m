(* ::Package:: *)

(* ::Title::Closed:: *)
(*core*)


(* ::Section::Closed:: *)
(*Dynamic Mode Decomposition*)


(* ::Subsubsection::Closed:: *)
(*Pruning a list of singular values*)


(* ::Item:: *)
(*Given sigmas, *)


(* ::Subitem:: *)
(*Use the tolerances to discard small singular values.*)


(* ::Subitem:: *)
(*Choose the first r  values when no tolerances are given.*)


sigthresholder[sigmas_,{sigtol_,abstol_}]:=Module[{threshold},
(* Relative criterion ONLY *)
threshold = sigmas[[1]]*sigtol;
 Select[sigmas,(#>= threshold)&]
];


sigthresholder[sigmas_,r_]:=Take[sigmas,r];


(* ::Subsubsection::Closed:: *)
(*Cut down U,V to match the reduced \[CapitalSigma]*)


thinsvd[usv_,sigmas_]:=Module[{nsigmas = Length@sigmas},
{Take[(usv[[1]])\[ConjugateTranspose],nsigmas]\[ConjugateTranspose],DiagonalMatrix[sigmas],Take[(usv[[3]])\[ConjugateTranspose],nsigmas]\[ConjugateTranspose]}
];


(* ::Subsubsection::Closed:: *)
(*Condition number and tail lost in reduced SVD*)


getcondandtail[sigmas_,n_]:=Module[{len =Length@sigmas},{If[sigmas[[-1]] > $MinMachineNumber,sigmas[[1]]/sigmas[[-1]],sigmas[[1]]/$MinMachineNumber],If[len > n,sigmas[[n+1]],0]}];


(* ::Subsubsection::Closed:: *)
(*Thin SVD + Performance metrics*)


(* ::Item:: *)
(*Find thin SVD of x*)


(* ::Subitem:: *)
(*Threshold singular values with a relative criterion*)


(* ::Item:: *)
(*Include the condition number and largest singular value dropped in thresholding*)


pickysvd[x_,sigtols_]:=Module[{usv,threshold,rawsigmas,sigmas,nsigmas,condandterror},
(* Get SVD of x *)
(* 
UpTo : To compute all Subscript[\[Sigma], i] so that we can threshold ourselves
*)
usv = SingularValueDecomposition[x,UpTo@Min[Dimensions[x]] ];
(* Here we choose the singular values that are considered and those that are neglected *)
(* Essentially, in addition to the relative tolerance implemented via sigtol, you also put an absolute threshold *)
(* I think the absolute threshold was put because it was ignoring a lot of zeros ???? I don't remember very well *)
rawsigmas = Diagonal[usv[[2]]];
sigmas = sigthresholder[rawsigmas, sigtols];
(* Find unthresholded \[Kappa](x) and the largest Subscript[\[Sigma], i] not included *)
condandterror = getcondandtail[rawsigmas,Length@sigmas];
(* Choose the corresponding rows of U^H and V^H *)
{thinsvd[usv,sigmas],condandterror}
];



(* ::Subsubsection::Closed:: *)
(*Exact DMD*)


coredmd[x_,y_,addons_,sigtols_,flavour_]:=Module[{usv,condandterror,apsuedo,redcond,cte},

{usv ,condandterror}= pickysvd[N@x,sigtols];
apsuedo = Dot[N@ y,usv[[-1]], (Inverse[usv[[2]]]),usv[[1]]\[ConjugateTranspose]];
(* Include the condition number of the reduced xmatrix*)
redcond = usv[[2,1,1]]*carefulinverse[usv[[2,-1,-1]]];
cte = Riffle[condandterror,redcond];
{apsuedo,cte}

]/;flavour=="direct"


(* ::Subsubsection::Closed:: *)
(*Companion DMD : Barebones*)


(* ::Item:: *)
(*Optimal linear 1 - step predictor & some error info*)


naivegetcstar[data_,sigtols_]:=Module[{xmat,xfinn,n,xdagger,cte},
xmat =( Most[data\[Transpose]])\[Transpose];
xfinn = data[[All,-1]];
n = Length@((xmat)\[Transpose]);
{xdagger,cte} = coredmd[xmat,IdentityMatrix[n],{},sigtols,"direct"];
{xdagger.xfinn,cte}
];



(* ::Section::Closed:: *)
(*Lifting the observations*)


(* ::Text:: *)
(*Not well documented, although tested.*)
(*Use of coordinates as observables obviates much of this infrastructure*)
(**)


(* ::Subsection::Closed:: *)
(*Dictionaries*)


(* ::Subsubsection::Closed:: *)
(*Polynomial *)


getflist[dindex_,maxpolydegree_,x_,basishead_]:=Table[Derivative[0,dindex][basishead][i,x],{i,0,maxpolydegree}];


varsep2whole[bases1d_]:=((Flatten@Outer[Times,##])&)@@bases1d;


psi[x_,{mpdegree_,basishead_}]:=
(* Rest ensures that the constant function is NOT used.*)Rest[varsep2whole@(Map[getflist[0,mpdegree,#,basishead]&,x])];



(* ::Subsubsection::Closed:: *)
(*Thin plate spline*)


splinevals[x_,basepts_]:= ((#^2*Log[#])&)[Map[Norm,((basepts\[Transpose] - x)\[Transpose])]];


(* ::Subsection::Closed:: *)
(*Generic lifting + makeextended*)


needpolytri[flavour_]:=Switch[flavour,
"polytri",True,
"delay",True,
_,False];


insideQ1D[x_,{a_,b_}]:=OrderedQ[{a,x,b}]&&x!= b;
insideQnD[x_,boxranges_]:=And@@MapThread[(insideQ1D[#1,#2])&,{x,boxranges}];



boxevals[x_,{widths_,numvec_,nboxes_,boxvol_,boxranges_,trueboxranges_}]:=If[insideQnD[x,trueboxranges],naiveboxevals[x,{widths,numvec,nboxes,boxvol,boxranges,trueboxranges}],SparseArray[{},{nboxes}]];


(* ::Item:: *)
(*Wrapper for psi and splinevals*)


lift[x_,params_,flavour_]:=psi[x,params]/;flavour=="polybox"
lift[x_,params_,flavour_]:=splinevals[x,params]/;flavour=="spline"
lift[x_,params_,flavour_]:=boxevals[x,params]/;flavour=="boxfuns"


(* ::Subsubsection::Closed:: *)
(*Monomials..(efficiently ?)*)


(* ::Item:: *)
(*Appends the desired monomial lifts to mat*)


(* ::Item:: *)
(*If no powerlist, returns the input*)


pascalift[mat_,{powlist_}]:=If[(* Do you want something beyond linear functions ? *)Length[powlist]!= 0 ,(* Append them to the given matrix *)Join[mat,Map[Times@@MapThread[(splpower[#1,#2])&,{mat,#}]&,powlist]],(* Else, keep what you were given *)mat];


(* ::Subsubsection::Closed:: *)
(*User-specified functions :*)


(* ::Item:: *)
(*Must behave like psi, splinevals*)


(* ::Item:: *)
(*Eg: Constant function*)


userlift[mat_,uparams_,userhead_]:=If[Length[userhead]!=0,Map[userhead[#,uparams]&,mat\[Transpose]]\[Transpose],{}];


(* ::Text:: *)
(*Lifts each delayed coordinate separately*)


superbasislift[mat_,params_,flavour_,obsdim_]:=Map[basislift[#,params,flavour]&,Partition[mat,obsdim]];


getXtendedmats[mat_,params_,flavour_,uparams_,userhead_,obsdim_]:=Module[{basisliftedrows,uliftedrows},
basisliftedrows = superbasislift[mat,params,flavour,obsdim];
uliftedrows = Map[userlift[#,uparams,userhead]&,Partition[mat,obsdim][[{1,2}]]];
MapThread[Flatten[Prepend[#1@basisliftedrows,#2],1]&,{{Most,Rest},uliftedrows}]
];


(* ::Subsubsection::Closed:: *)
(*splscale*)


splscale[scales_,mat_]:=
(* Filter out the indices where the coordinate gets multiplied by 0 *)
(* This eliminates 0 rows in the data matrix *)
Times@@Map[Pick[#,Clip[scales],1]&,{scales,mat}];


(* ::Subsubsection::Closed:: *)
(*Creates detailed I/P from simple I/P within paramsplus*)


trigscalelts2list[nlist_]:=Module[{nmax,m0,rawlist,checklist},
nmax = Max@nlist;
m0 = Length@nlist;
(* Take all possible combinations of all scales you need *)
rawlist = Flatten[(Outer[List,##]&)@@(Map[Range[0,#]&,nlist]),m0-1];
(* Filter out the ones that violate the max permissible scale sum :  ntotal *)
checklist = Map[TrueQ[Total[#]<= nmax]&,rawlist];
Rest@Pick[rawlist,checklist]
];
