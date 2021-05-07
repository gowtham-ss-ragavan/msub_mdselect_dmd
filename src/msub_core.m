(* ::Package:: *)

(* ::Title:: *)
(*core*)


(* ::Subsubsection::Closed:: *)
(*Distance to spectral equivalence with DFT : \[Delta]^\[Mu]_{rel}*)


obs4c[c_]:=Norm[c-(-1)]/Sqrt[Length[c]];


(* ::Subsubsection:: *)
(*Plot Subscript[\[Rho], DFT]for a range of delays + RoM order*)


(* ::Subsection:: *)
(*Companion matrix*)


(* ::Subsubsection:: *)
(*Generation*)


getCompanionmat[c_]:=Transpose@Rest@Append[IdentityMatrix[Length[c]],c];


(* ::Subsubsection:: *)
(*Spectrum *)


getCompanionEvals[c_]:=Eigenvalues[getCompanionmat[c]];


(* ::Subsubsection:: *)
(*Eigen-decomposition *)


getCompanionEigsys[c_]:=Module[{localn = Length[c],vals,scales,vecs},
vals = getCompanionEvals[c];
vecs =getCompanionEvecs[vals,localn];
{vals,vecs}
];


(* ::Subsubsection:: *)
(*Coefficients of the monic polynomial with given roots*)


roots2coeffs[roots_]:=Module[{x},CoefficientList[Times@@(x-roots),x]];


(* ::Subsection:: *)
(*Vandermonde & related operations*)


(* ::Subsubsection:: *)
(*Vandermonde matrix from specified nodes*)


(* ::Item:: *)
(*n: The highest power to raise the nodes to *)


(* ::Item:: *)
(*\[CapitalTheta] = getVandermonde[{Subscript[\[Lambda], 1],Subscript[\[Lambda], 2],...Subscript[\[Lambda], r]},n] *)


(* ::Item:: *)
(*O/P ~  #(nodes) x (n+1) *)


getVandermonde[nodes_,n_]:= Table[nodes^i,{i,0,n}]\[Transpose] ;


(* ::Subsubsection::Closed:: *)
(*Delayed Vandermonde matrices*)


getZmat[nodes_,n_,ndelays_]:=Join@@NestList[(DiagonalMatrix[nodes].#)&,getVandermonde[nodes,n],ndelays];


(* ::Subsubsection::Closed:: *)
(*Subtracting a "generalized" mean*)


(* ::Item:: *)
(*If \[Lambda] = 0, return the data set as is*)


(* Compute the mean subtracted data matrix *)
getmusubeddatamat[data_,\[Lambda]_]:=Module[{n,gpvec},
n = Length[data\[Transpose]]-1;
gpvec = getgpvec[\[Lambda],n];
data - {(data.(gpvec^-1))*(1/(n+1))}\[Transpose].{gpvec}
];


(* ::Subsubsection::Closed:: *)
(*Remove all elements that are close to a specified number*)


abstractremovelambfromset[set_,setpicker_,lamb_,tolerance_]:=Pick[set,(* Elements that are far eoungh will be True *)Thread[Abs[setpicker[set]-lamb]>tolerance]];


removelambfromset[set_,lamb_,tolerance_]:=abstractremovelambfromset[set,#&,lamb,tolerance];


(* ::Section:: *)
(*Constructing le PDFs*)


(* ::Subsubsection::Closed:: *)
(*Un-delaying a DMD mode*)


(* ::Item:: *)
(*Removes the effect of taking delays from a DMD mode*)


(* ::Subitem:: *)
(*Computationally, this is GLA wrt the corresponding DMD eigenvalue*)


(* ::Subitem:: *)
(*Helps when you want to compare DMD modes across a range of delays*)


(* ::Subsubitem:: *)
(*Justified only when you do no mean-subtraction*)


delayfilter[rawmode_,eval_,crows_,ndelays_]:=Module[{modemat,rawvec,hitvec},
modemat = Partition[rawmode,crows]\[Transpose]; (* crows x ndelays+1*)
rawvec = getgpvec[carefulinverse[eval],ndelays];
hitvec = rawvec/(ndelays+1);
modemat.hitvec
];


(* ::Subsubsection::Closed:: *)
(*Rating DMD eigenvalues*)


(* ::Item:: *)
(*Computes a score* describing the relative weight of each DMD eigenvalue*)


(* ::Subitem:: *)
(*Currently,  Mode norm/ Sum of mode norms*)


(* Has been tested for when jac was actually used in computing probs *)
cmodesvals2pdf[c_,modes_,vals_]:=Module[{wts,jacs,probs,npow = 2,dpow = 1},
wts = Map[Norm,modes];
(*jacs = clamb2jacs[c,vals];*)
probs = ((#/Total[#])&)@(wts^2)(*wtsnjacs2rfactors[wts,jacs,{npow,dpow}]*);
{vals,probs}
];


(* ::Subsubsection::Closed:: *)
(*Data-scored eigenvalues of a Companion matrix*)


(* ::Item:: *)
(*This part...can go awfully awry when there are rogue eigenvalues*)


cX2valswts[c_,xmat_,filterchoice_,crows_,ndelays_]:=Module[{vals,vecs,rawmodes,modes},
(* Eigen-decomposition of Companion matrix *)
{vals,vecs}=getCompanionEigsys[c];
(* Find the DMD-modes *)
rawmodes = (xmat.(vecs\[Transpose]))\[Transpose];
(* Un-delay the modes, if filter is a go *)
modes= If[filterchoice=="yes", MapThread[delayfilter[#2,#1,crows,ndelays]&,{vals,rawmodes}],rawmodes];
cmodesvals2pdf[c,modes,vals]
];


(* ::Subsubsection::Closed:: *)
(*Optimal inclusive polynomial*)


(* ::Text:: *)
(*Find the polynomial whose coefficients *)


(* ::Subitem:: *)
(*Contain  desired eigenvalues : targetevals*)


(* ::Subitem:: *)
(*Closest to a given polynomial: cstar*)


closestructvec[cstar_,targetevals_,sigtols_]:=Module[{pointedc0,cstarlen,c0len,cstaronboots,locrules,locxmat,locyvec,locdata,alpha,cte,cstarmpointed,spuriousvals,structedcvec,loczmat,evalstar,subsetnum},
(* Compute the min bubble size to draw about eigenvalues of CompanionMatrix from cstar, to contain targetevals *)
evalstar = getCompanionEvals[cstar];
subsetnum=pert4lhsinrhs[targetevals,evalstar];
(**)
pointedc0 = roots2coeffs[targetevals];
{cstarlen,c0len} = Map[Length,{cstar,targetevals}];
If[c0len >= cstarlen (* Even if you have as many DoFs (cstarlen) as the number of constraints (cstarlen), there is nothing to be done because we need to keep the coefficient monic as in the higher degree case *),
{subsetnum,-Most[pointedc0],{}},
loczmat = ToeplitzMatrix[PadRight[pointedc0,cstarlen],PadRight[{pointedc0[[1]]},cstarlen-c0len+1]];
locxmat = loczmat[[All,1;;-2]];
locyvec = cstar + PadLeft[Most[pointedc0],cstarlen];
locdata = Append[locxmat\[Transpose],locyvec]\[Transpose];
{alpha,cte} = naivegetcstar[locdata,sigtols];
structedcvec = loczmat.Append[alpha,-1];
spuriousvals = (getCompanionEigsys[alpha])[[1]];
{subsetnum,structedcvec,spuriousvals}
]
];


(* ::Subsubsection::Closed:: *)
(*DMD-eigenvalues expected from theory*)


(* ::Item:: *)
(*Implements Table 6.1 in draft_jan20 *)


getevalpartitions[truevals_,nfunda_,n_,rate2sub_,restols_]:=Module[{truesansrate,truewithrate,abvanilla},
truesansrate = removelambfromset[truevals,rate2sub,restols[[1]]];
truewithrate = Append[truesansrate(* Not truevals to avoid a potential double eigenvalue in structedcpres *),rate2sub];
abvanilla=
Switch[Mod[n+1,nfunda+1],
0,(* Case 1 *) 
Insert[(*ms and msdel won't have \[Lambda] in \[Sigma](T[c]) as it is a clean cut in Case 1*)ConstantArray[{truesansrate,{}},2],{truesansrate,{rate2sub}} (* mspres has \[Lambda] in \[Sigma](T[c]) as it is explicitly preserved, but gets knocked out when constructing the PDF *),2],
_,(* Case 2 *)
ConstantArray[(* \[Lambda] is always present in Subscript[\[Sigma], num] for Case 2*){truewithrate,{}},3]
];
Join[(* Vanilla remains the same for both cases: We know that the true eigenvalues must be within the Subscript[\[Sigma], num]*){{truevals,{}}},abvanilla]
];


(* ::Subsection:: *)
(*Additional abstractions within getcstar++ :*)


(* ::Subsubsection::Closed:: *)
(*Are the leading few elements > the rest ?*)


(* ::Item:: *)
(*Distinguishes the two cases for \Gamma[Z,B,B_0,c]*)


non0coresQ[probs_,ncores_]:=Module[{cores,noncores},
{cores,noncores} = TakeDrop[probs,ncores];
TrueQ[Min[cores] > Max[noncores]]
];


(* ::Subsubsection:: *)
(*Numerical checks for mode-selection theorems*)


(* ::Item:: *)
(*Gives {\[Rho]_{Subset},\[Delta]_{Trivial}}*)


numchks4msubthrms[xmat_,evalparts_,cvec_,sigtols_]:=Module[{n,corevals,auxevals,structedc,spuriousvals,ncorevals,structevals,structevecs,structedpdf,rperror,respercent},
(*n is locally computed because updatednumchks handles all 4 cases wherein n differs for the sub-then-delay case *)
n = Length@cvec;
{corevals (*B*),auxevals} = evalparts;
(* Find the closest c that contains all the specified eigenvalues in evalparts *)
{rperror (*Subscript[\[Rho], subset]*),structedc(* ONLY WHEN YOU HAVE SPACE TO ACCOMODATE ALL YER COSNTRAINTS !!!!!*),spuriousvals} = closestructvec[cvec,(*Subscript[B, 0]*)Join[corevals,auxevals],sigtols];
ncorevals = Length@corevals;
structevals = Join[corevals,auxevals,spuriousvals];
respercent = If[Length[structevals]>n,
(* #[Constraints] > #[DMD order] \[Equal]> Cannot construct a PDF *)1,
(* Else, strcutedc will have just the right dimensions to hit xmat from the right*)structevecs = getCompanionEvecs[structevals];
structedpdf = cmodesvals2pdf[structedc ,(xmat.(structevecs\[Transpose]))\[Transpose],structevals];
(*---------- Compute the unexplained fraction --------------------------------*)
1 -  Total[Take[structedpdf[[2]],ncorevals]]
];
{rperror,respercent}
];


(* ::Subsubsection::Closed:: *)
(*Prepare for DMD & mean-subtracted cousins*)


(* ::Item:: *)
(*Process the lifted time-series*)


(* ::Subitem:: *)
(*Remove the mean*)


(* ::Subitem:: *)
(*Delay if needed*)


(* ::Item:: *)
(*Choose to filter modes for vanilla DMD*)


(* ::Item:: *)
(*Generate the set of eigenvalues you expect to be captured in each flavour*)


prep4cases[data_,\[Lambda]_,truevals_,restols_,nfunda_]:=Module[{zmat,zlambdams,zlambdamsdel,n,evalparts,flavours,dataorgs},
zmat =  data;
zlambdams = getmusubeddatamat[zmat,\[Lambda]];
zlambdamsdel = Join@@Map[#[zlambdams\[Transpose]]\[Transpose]&,{Most,Rest}];
n = Length[zmat\[Transpose]]-1;
evalparts = getevalpartitions[truevals,nfunda,n,\[Lambda],restols];
flavours = {"vanilla","ms","mspres","msdelay"};
dataorgs = {zmat (* Usual *),zlambdams (* Mean subtracted *),zlambdams (* Remove mean - Preserve \[Lambda]*),zlambdamsdel(* Remove mean, then delay *)};
{dataorgs,evalparts,flavours}
];


(* ::Subsubsection::Closed:: *)
(*Companion DMD & mean-subtracted versions *)


(* ::Item:: *)
(*In addition to usual regression results, it gives*)


(* ::Subitem:: *)
(*A list of scored DMD-eigenvalues*)


(* ::Subitem:: *)
(**)


(* ::Item:: *)
(*Works best for systems of dimensions > 2*)


getcstar[data_,\[Lambda]_,sigtols_,restols_,evalparts_,flavour_,crows_,ndelays_]:=Module[{xmat,m,n,cvec,cte,pdf,ratereps,bandmat,rawcvec,corevals,auxevals,spuriousvals,rperror,structedc,ncorevals,structevals,structedpdf,structevecs,respercent,filterchoice}, 
xmat =zmat2xmat[data];
{m,n} = Dimensions[xmat];
If[flavour != "mspres",
(* Simple regression if you are not preserving \[Lambda] *)
{cvec,cte} = naivegetcstar[data,sigtols],
(* Encode \[Lambda] preservation by searching for c in an appropriate subspace *)
bandmat = SparseArray[{Band[{1,1}]-> -\[Lambda],Band[{2,1}]-> 1},{n+1,n}];
{rawcvec,cte} = naivegetcstar[data.bandmat,sigtols];
cvec =Most[bandmat.Flatten[{rawcvec,-1}]]
];
(* try to recover the true modes when you are vanilla, in which case you know the delayed mode is made of iterative scaling of the true mode *)
filterchoice = If[flavour=="vanilla","yes","no"];
If[(* 3D and above are cool : Case 1 2D systems (not minn) will not have respercent and rperror computed as they are indistinguishable from No info systems *)Length[Flatten[evalparts]] <= 1,(* If nothing is known about the spectrum, simply compute the PDF using cvec, and set rperror,respercent to HIGH values  *)
ratereps = {cX2valswts[cvec,xmat,filterchoice ,crows,ndelays],1,1},
{rperror,respercent} = numchks4msubthrms[xmat,evalparts,cvec,sigtols];
ratereps=  {{},rperror,respercent};
If[flavour=="vanilla",ratereps[[1]] = cX2valswts[cvec,xmat,filterchoice ,crows,ndelays]]
];
{cvec,ratereps,cte}
];


(* ::Subsubsection:: *)
(*Analyse a trajectory with DMD & mean-subtracted versions*)


getcquad[data_,\[Lambda]_,sigtols_,restols_,truevals_,nfunda_,crows_,ndelays_]:=Module[{datapartsflavs},
(* First create the various inputs for each flavour of Companion DMD *)
datapartsflavs = prep4cases[data,\[Lambda],truevals,restols,nfunda];
MapThread[getcstar[#1,\[Lambda],sigtols,restols,#2,#3,crows,ndelays]&,datapartsflavs]\[Transpose](* [[Attributes, Cases ]]*)
];


(* ::Section:: *)
(*Sneaky estimation of rootoptimalsubsets*)


(* ::Subsubsection:: *)
(*Score c* using its error diagnostics*)


(* ::Item:: *)
(*Lower redcond and te are favoured*)


(* ::Subitem:: *)
(*redcond : Condition number of the low rank approximation of X (<=  a user specified value)*)


(* ::Subitem:: *)
(*te : Largest singular omitted in the above reduction*)


cte2wt[{cond_,redcond_,te_},{cpow_,tepow_}]:=1/(1+redcond^cpow*te^tepow);


(* ::Subsubsection:: *)
(*Remove trailing zeros in a vector*)


maketrim[vec_,restols_]:=Module[{n,ttable,lastindex},
n = Length@vec;
ttable = Map[TrueQ[#>restols[[2]]]&,Abs@vec];
lastindex = (Pick[Range[n],ttable])[[-1]];
Take[vec,lastindex]
];


(* ::Subsubsection:: *)
(*Roots of the polynomial represented by a vector*)


(* ::Item:: *)
(*Tries to find the one of the lowest degree*)


(* ::Text:: *)
(*!!!CANNOT HANDLE THE ZERO VECTOR !!!*)


coeffs2roots[coeffs_,restols_]:=Module[{trimcoeffs},
trimcoeffs = maketrim[coeffs,restols];
Switch[Length[trimcoeffs],
1,Nothing,
_,getCompanionEvals[Most[trimcoeffs]/-trimcoeffs[[-1]]]
]
];


(* ::Subsubsection:: *)
(*Evaluating a quadratic form on GPs from a discrete set*)


(* ::Item:: *)
(*Used to find Subscript[\[Rho], SoS][\[Lambda]]^2 for \[Lambda] \[Epsilon] roots*)


graderootcandidates[dmat_,roots_,restols_,(* Length@dmat *)n_]:=Module[{vecs,qforms,restests},
vecs = Map[getgpvec[#,n-1]&,roots];
qforms = Map[(* Technically always non-ngeative...but we have some pesky small imaginary parts due to the numerics...*)Abs@(*Subscript[\[Rho], SoS][\[Lambda]]^2 = Subscript[gp, n+1][\[Lambda]]^H dmat Subscript[gp, n+1][\[Lambda]]*)(Dot[Conjugate[#],dmat,#])&,vecs];
{roots,qforms}
];


(* ::Subsubsection:: *)
(*DMD-eigenvalues from a weighted list of c**)


(* ::Item:: *)
(*~ Steps 3-5 in Algorithm 6.1*)


aw2sneakyropts[aw_,restols_,sigtols_]:=Module[{dmat,localn,usv,cte,numericu2,umat,u2,pdf,roots,probs,qforms},
(* dmat = (Overscript[\[Alpha], -] Overscript[W, -])(\[Alpha] W)^T = Conjugate[aw].aw\[Transpose] *)
dmat = Total@Map[KroneckerProduct[Conjugate[#],#]&,aw\[Transpose]];
localn = Length@dmat;
(* Get Subscript[U, 2] : A basis for Ker[aw\[ConjugateTranspose]]*)
{usv,cte}=pickysvd[aw,sigtols];
(* Find the roots for this particular polynomial *)
roots = coeffs2roots[(* Get the largest left singular vector *)usv[[1,All,1]],restols];
(* Grade it : Computing Subscript[\[Rho], SoS][\[Lambda]]^2 for roots *)
qforms = (graderootcandidates[dmat,roots,restols,localn])[[2]];
{roots,qforms}
];


(* ::Subsubsection:: *)
(*Sort a list of pairs according to the second coordinate*)


orderqfpairs[rqfpairs_]:=Module[{orderpls},
(* Get Ordering for the second column *)
orderpls =Ordering[rqfpairs[[All,2]]];
(* Use to rearrange the whole *)
rqfpairs[[orderpls]]
];


(* ::Subsubsection:: *)
(*DMD eigenvalues from a list of c**)


(* ::Item:: *)
(*Implements  Algorithm 6.1*)


(* ::Subitem:: *)
(*rmin: "r"*)


(* ::Subitem:: *)
(*chosendeg:"n"*)


(* ::Item:: *)
(*Generates weights, calls aw2sneakyropts and then gives you its take on the estimates*)


sneakyestimatetruevals[vals_,crunchgrub_,opstyle_]:=sneakyestimatetruevals[vals,crunchgrub[chosendeg],crunchgrub[testdegs],{crunchgrub[sigtols],crunchgrub[restols]},opstyle,crunchgrub[cpow],crunchgrub[tepow],crunchgrub[rmin]];


sneakyestimatetruevals[vals_,chosendeg_,testdegs_,{sigtols_,restols_},opstyle_,cpow_,tepow_,rmin_]:=Module[{chosenindex,localwts,wtsum,localcmat,localalpha,localphaW,sneakyropts,rqfpairs,estruevals,objvalplots,specplot,chosenspecplot},
(* Convert chosendeg to its index wrt testdegs*)
chosenindex = ( Position[testdegs,chosendeg])[[1,1]];
(*---------------- Assemble weight matrix W -----------------------*)
(* Compute the weights associated with each column of \[Alpha] by using the condition number and tail error in computing c via coredmd *)
localwts = Map[cte2wt[#,{cpow,tepow}]&,vals2vanillactes[vals,chosenindex]];
wtsum = Total[localwts];
localwts*= carefulinverse[wtsum];
(*---------------- Assemble cvec matrix \[Alpha] -----------------------*)
(* Add -1 below cvecs to make \[Alpha] *)
localcmat = (vals2vanillacvecs[vals,chosenindex])\[Transpose];
localalpha = dresscmat[localcmat];
(*---------------- Assemble \[Alpha] W -----------------------*)
(* Weigh each column of \[Alpha] and use aw2sneakyropts *)
localphaW = (localalpha\[Transpose]*localwts)\[Transpose];
(*---------------- Get eigenvalue estimates: Overscript[\[Nu], ^] -----------------------*)
sneakyropts (*{roots,qforms}*) = aw2sneakyropts[localphaW,restols,sigtols];
(* Convert into a list of root-QF(Subscript[\[Rho], SoS][root]^2) pairs *)
rqfpairs = (* Sort in the order of increasing Subscript[\[Rho], SoS][root]^2 *)orderqfpairs[(* root-QF pairs *)sneakyropts\[Transpose]];
(* Pick the best rmin pairs *)
estruevals = rqfpairs[[ Range[Min[rmin,Length[rqfpairs]]](* In case the singular vector decides to prune the roots beforehand ...*),1]];
(*O/P as desired*)
Switch[opstyle,
"brevity",estruevals,
_,(* Plot the variation of Subscript[\[Rho], SoS][\[Lambda]]^2 with \[Lambda] *)
objvalplots =ListPlot[splog10[rqfpairs[[All,2]]],PlotRange->All,AxesLabel-> {"\[Lambda] index","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(\[Rho]\), \(SoS\)]\)[\[Lambda]\!\(\*SuperscriptBox[\(]\), \(2\)]\)]"},ImageSize->Large]; 
(* Different visualization fot he above for rqfpairs : Subscript[\[Rho], SoS]^2 is represented with colour *)
specplot = powplaneplot[rqfpairs,(-splog10[#])&];
(* Plot estruevals to see where they eventually fall *)
chosenspecplot = ListPlot[complexto2D@estruevals,PlotRange->All,PlotStyle->Directive[Black,PointSize[Large]],PlotMarkers->"O"];
Print[objvalplots];
Print[Show[specplot,chosenspecplot]];
Print["\!\(\*SubscriptBox[\(\[Rho]\), \(SoS\)]\)-- as Red \[LongRightArrow] Blue"];
{rqfpairs,estruevals}
]
];
