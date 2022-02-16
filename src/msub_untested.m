(* ::Package:: *)

(* ::Title:: *)
(*untested*)


(* ::Subsubsection::Closed:: *)
(*Explicit computation of right eigenvectors of a Companion matrix*)


getCompanionEvecs[vals_,clen_]:=Module[{vecs,scales},
vecs = Table[roots2coeffs[Drop[vals,{i}]],{i,clen}];
scales = MapThread[(getgpvec[#1,clen-1] . #2)&,{vals,vecs}];
(scales^(-1))*vecs
];


getCompanionEvecs[vals_]:= getCompanionEvecs[vals, Length[vals]];


(* ::Subsubsection::Closed:: *)
(*Update the numerical checks for mode-selection theorems using estimatruevals*)


(* ::Item:: *)
(*When the true Koopman eigenvalues are unknown, one needs to estimate them before computing the numerical checks*)


updatednumchks[data_,\[Lambda]_,truevals_,restols_,sigtols_,nfunda_,locquads_]:=Module[{dataorgs,evalparts,flavours},
{dataorgs,evalparts,flavours} = prep4cases[data,\[Lambda],truevals,restols,nfunda];
MapThread[numchks4msubthrms[zmat2xmat[#1],#2,#3,sigtols]&,{dataorgs,evalparts,locquads}]
];


(* ::Subsubsection:: *)
(*Drop unorganized elements of the list of associations named vals*)


onlythequads[vals_] := vals[[All,(* Lose the other Keys as they may not ahve as man elements as in these three *){Key[cquads], Key[ratequads], Key[equads]}]];


(* ::Subsubsection:: *)
(*Generating data to help choose \hat{r} : The approximate dimension of the underlying invariant subspace*)


(* ::Subitem:: *)
(*Compute the pair-wise \rho_{Subset} over the range of "n" under purview*)


(* ::Subitem:: *)
(*Will then need to UpperTriangularize, as we are interested in the smaller being contained in the larger set*)


getetvdepwrtcdeg[vals_, cgrubin_] := Module[{cgrub = cgrubin, localetvs, pairdts, plotcoords, zcoords, locvals},
	(*Common eigenvalues,over ICs and permitted delays,at varying choices of n*)
	localetvs = 
	Table[(
		AssociateTo[cgrub, {chosendeg -> i,(*Want all the eigenvalues for a chosen value of n (chosendeg)*)rmin -> i, delaymin -> i - 1}];
        locvals = keepthemgoodelays[onlythequads[vals], cgrub]; 
        sneakyestimatetruevals[locvals, cgrub, "brevity"]
        ), 
   {i, cgrubin[testdegs]}];
   (*Want to know how much the low dimensional estimates are contained in the high dimensional versions*)(*So,first compute the pair-wise values of \rho_{Subset}*)
   pairdts = Outer[pert4lhsinrhs, localetvs, localetvs, 1];
   plotcoords = Flatten[Outer[{##} &, cgrubin[testdegs], cgrubin[testdegs]], 1];
   (*Assuming the first coordinate represents the first argument,we zero out all those with a larger first argument,AFTER taking a log10*)
   zcoords = Flatten@UpperTriangularize@splog10[pairdts];
   ldplot[plotcoords, zcoords, {"Bigger n", "Smaller n"}, {(*Skip the diagonals*)-12, -0.01 (*Skip the zeros*)}, "TemperatureMap"(* A more vivid color scheme than the default*)]
   ];


(* ::Section:: *)
(*Mean sub mods*)


(* ::Subsubsection::Closed:: *)
(*Lift points iff they are within an origin centered hypercube*)


clippedbasislift[mat_,params_,flavour_]:=Module[{logiclist,insideones,liftedones,outlifted,index},
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


(* ::Subsubsection::Closed:: *)
(*MOD: Lift indiscriminately*)


(* ::Item:: *)
(*naivebasislift for now*)


(* ::Subitem:: *)
(*Change to clippedbasislift if doing S-DMD*)


(* ::Item:: *)
(*params[[-1]]: \tilde{C} ~ Coefficients of the dictionary in terms of the eigenfunctions of the underlying Koopman invariant SS*)


basislift[mat_,params_,flavour_]:=params[[-1]] . (naivebasislift[mat,Most[params],flavour]);


(* ::Subsubsection::Closed:: *)
(*Parameter generation for different choices of dictionaries*)


(* ::Item:: *)
(*msub specific part: Generation of a random coefficient matrix \tilde{C}*)


getparamsplus[dummyvars_,flavour_,flavourparams_,boxflavourparams_,obsdim_,cmatseed_,crows_,crank_]:=Module[{params,m15,cmat},
params ={ Switch[flavour,
"delay",{{}},
"spline",getrandpoints[flavourparams],
"polybox",flavourparams,
"polytri",{Flatten[Table[Keys@CoefficientRules[(Plus@@dummyvars)^i,dummyvars],{i,2,flavourparams (* Erstwhile maxpolydegree *)}],1]},
"trigtri",{trigscalelts2list[flavourparams]},
"boxfuns",cookboxgrub[flavourparams(*{boxranges,centercounts}*)],
_,Print["Exiting as flavour NOT specified !!!!"];Abort[]
],(cookboxgrub[boxflavourparams])[[-1]] (* Used to select which points to exclude before the lift *)};
m15 = Length@ (( (* This will be unevaluated as 1 is not a matrix...So, we separate them and use only the last element...*)List@@basislift[{ConstantArray[0,obsdim]}\[Transpose],Append[params,1],flavour])[[-1]]);
cmat = If[cmatseed == {}(* If unspecified, create using getrandmat *), usvdot@((* randmat's OP is in usv form*)getrandmat[crows,m15,{crank,{}}]),cmatseed (* Use given Cmatrix*)];
(* Include as final parameter in params *)
{AppendTo[params,cmat],m15}
];




getparamsplus[dummyvars_,obsdim_,liftgrub_]:=getparamsplus[dummyvars,liftgrub[flavour],liftgrub[flavourparams],liftgrub[boxflavourparams],obsdim,liftgrub[cmatseed],liftgrub[crows],liftgrub[crank]];


(* ::Subsection:: *)
(*meansuboneshot*)


(* ::Subsubsection:: *)
(*No structing*)


(* ::Item:: *)
(*The "main" function for performing all computation relevant to generating the numerical checks *)


(* ::Item:: *)
(*Prepare - Lift - DMD - Organize*)


meansuboneshot[ndelays_,tinit_,vfield_,rawics_,chosenputty_,sigtols_,restols_,flavour_,flavourparams_,construct_,basistruct_,userparams_,rawuserlist_,invsamplesets_,ipfuns_,opfuns_,testinvsamplesets_,testipfuns_,testopfuns_,boxflavourparams_,cmatseed_,crows_,crank_,testdelays_,testdegs_,rate2sub_,nfunda_,truevals_,ipgrub_]:=Module[{sampackage,tfinn,ics,obsdim,dummyvars,params,userhead,m1,m2,commoninvgrub,giveninvgrub,e0,f0,q0,r0,teste0,testf0,testq0,testr0,testinvgrub,residue,kerrors,amat,conerrors,liftedX,liftedY,userlist,datamat,m15,tseries,cmat,opgrub,crequads,crefun,matcarvefun,deldegmat,newchks,newrqs},
(*-----------------------------*)
(*------PRELIMINARIES-------*)
(*-----------------------------*)
userlist = rawuserlist;
(* Create the ICs*)
tseries =rawics ;
(* Convert the time series into data matrices that can be used by getXtended *)
(* #[delays] effected is ndelays +1, so that it compensates for the loss of 1 in getXtended *)
datamat =tseries2bigdatamat[tseries,ndelays];
(* Parameters that we'll need later on*)
obsdim = Length[rawics[[1]]]; 
(* Number of observables for 0 delays: dummayvars +  Constraint generation when using a delay basis *)
(*-----------------------------*)
(*------FLAVOUR AND PARAMS-------*)
(*-----------------------------*)
dummyvars = Table[Subscript[x, i],{i,obsdim}];
{params,m15} = getparamsplus[dummyvars,flavour,flavourparams,boxflavourparams,obsdim,cmatseed,crows,crank];
(*-----------------------------*)
(*--------SWITCHES-------*)
(*-----------------------------*)
If[construct,AppendTo[userlist,(1)&]];
userhead = Through[userlist[##]]&;
m1 = Length@userlist; (* Number of user specified functions *) 
m2 = m15*(ndelays+1);
(* Number of basis functions taken : May need to evaluate basislift to get this*)
(*-----------------------------*)
(*-----GENERATE CONSTRAINT MATRICES-----*)
(*-----------------------------*)
commoninvgrub = {ndelays,dummyvars,params,flavour,userparams,userhead,construct,basistruct,m1,m2,obsdim,m15};
giveninvgrub = {invsamplesets,ipfuns,opfuns};
{liftedX,liftedY} =getXtendedmats[datamat,params,flavour,userparams,userhead,obsdim];
(* liftedY is of no use whatsover for this study since it simply uses Companion DMD*)
(*----------------------------------------*)
(* Compute c-r-e quads over the desired set of delays, trajectory length / system order *)
matcarvefun = (liftedX[[Range[m1 + (* The new face..NOT m15*)crows*(1+#1)],Range[(*(* Will do a delay, and so need atleast 3..and hence the 2*)2*)  (* testdegs are better written as testns... and so we need only take one additional column *)1+#2]]])&;
If[ipgrub==<||>(* Length[ipgrub] \[Equal] 0 *),
(*------------- First computational pass : No priors ----------------*)
crefun =getcquad[matcarvefun[#1,#2],rate2sub,sigtols,restols,truevals,nfunda,crows,#1]&;
crequads=Transpose[Outer[crefun,testdelays,testdegs],{2,3,1,4}];
opgrub=<|cquads-> crequads[[1]],ratequads-> crequads[[2]],equads-> crequads[[3]],mexparams->  {params,flavour,userparams,userhead,obsdim},i2bdmparams-> {vfield,tinit,tfinn,sampackage,chosenputty} |>,
(*------------- Updating numerical checks with new eigenvalue estimates ----------------*)
opgrub = ipgrub;
crefun =(* Compute rperror and respercent with truevals whiich is likely a product of estimatruevals*) updatednumchks[matcarvefun@@#1,rate2sub,truevals,restols,sigtols,nfunda,#2]&;
deldegmat = Outer[{##}&,testdelays,testdegs];
newchks = MapThread[crefun,{deldegmat,ipgrub[cquads]},2];
(*Update opgrub with newchks appropriately *)
newrqs = opgrub[ratequads];
newrqs[[(* Delays *)All,(* Degrees *)All,(* Cases *)All,{2,3}]]=newchks;
opgrub[ratequads]=newrqs
];
opgrub
];


(* ::Subsubsection::Closed:: *)
(*Structed*)


(* ::Text:: *)
(*UPDATE this list of keys !!!!*)


(* ::Item:: *)
(*Keys@trajgrub*)


(* ::Subitem:: *)
(*{discevals, tinit, tsamp, maxdelays, maxn, vfield, chosenputty, rawics}*)


(* ::Item:: *)
(*Keys@liftgrub*)


(* ::Subitem:: *)
(*{cmatseed, crank, crows, nprojs, flavour, flavourparams, boxflavourparams}*)


(* ::Item:: *)
(*Keys@priorsgrub*)


(* ::Subitem:: *)
(*{userparams, invsamplesets, ipfuns, opfuns, rawuserlist}*)


(* ::Item:: *)
(*Keys@testpriorsgrub*)


(* ::Subitem:: *)
(*{invsamplesets, ipfuns, opfuns}*)


(* ::Item:: *)
(*Keys@crunchgrub*)


(* ::Subitem:: *)
(*{savefile, testdelays, testdegs, sigtols, restols}*)


meansuboneshot[grubsNprior__Association]:=Module[{grubs,prior,args},
{prior,grubs} = TakeDrop[{grubsNprior},-1];
args = Join[{maxdelays,tinit,vfield,rawics,chosenputty,sigtols,restols,flavour,flavourparams,construct,basistruct,userparams,rawuserlist,invsamplesets,ipfuns,opfuns,testinvsamplesets,testipfuns,testopfuns,boxflavourparams,cmatseed,crows,crank,testdelays,testdegs,rate2sub,nfunda,truevals}/.(Normal[Join@@grubs]),prior];
meansuboneshot@@args
];



(* ::Section:: *)
(*Logistical spawn*)


(* ::Subsubsection::Closed:: *)
(*MOD: Random matrices with condition number = 1 when \[CapitalSigma] unspecified*)


getrandmat[m_, n_, {r_, sigmas_}] := Module[{buffer, sigs},
   buffer = SingularValueDecomposition[getfullrankmat[m, n]];
   sigs = Switch[Length@sigmas,
     0, ConstantArray[1, r],
     _, sigmas];
   thinsvd[buffer, sigs]
   ];


(* ::Subsubsection::Closed:: *)
(*Layering functions on 3D arrays*)


makecratesfunny[layer1fun_, layer2fun_, crates_] := ParallelMap[layer1fun, Transpose[Map[layer2fun, crates], {3, 1, 2}], {2}];


(* ::Subsubsection::Closed:: *)
(*Attach scalar attributes to a list of (x,y) coordinates*)


getxyzlist[xylist_, zlist_] := Append[xylist\[Transpose], zlist]\[Transpose];


(* ::Subsubsection:: *)
(*ListDensityPlot of coordinate list {(y,x)} with corresponding attributes {z}*)


(* ::Text:: *)
(*Generic version used to visualize variation of \rho_{Subset}[Estimate_{i},Estimate_{j}] as i,j range over the DMD model orders under study*)


ldplot[deldegcoords_,coordvalues_,framelabel_,plotrange_,colourfun_]:=Module[{xycoords,zcoords,xrange,yrange,regfun,ldpin},
xycoords=Reverse[deldegcoords,2];(*X:n,Y:#[delays]*)
zcoords=Flatten[coordvalues];
(*Generate the list of 3D coordinates requiredc for ListDensityPlot*)
ldpin=getxyzlist[xycoords,zcoords];
(* Sometimes ListDensityPlot extrapolates a bit outside what is given *)
(* Hence, we create a RegionFunction to do just that*)
{xrange,yrange} = Map[MinMax,xycoords\[Transpose]];
regfun= Function[{x,y,z},xrange[[1]] <= x <= xrange[[2]] && yrange[[1]] <= y <= yrange[[2]]];
(* Combine all options in the plot*)
ListDensityPlot[ldpin,InterpolationOrder->0,PlotLegends->Automatic,FrameLabel->framelabel,PlotRange->plotrange,PerformanceGoal->"Quality",LabelStyle->Directive[Black,15],RegionFunction->regfun,ColorFunction->colourfun]
];


(* ::Text:: *)
(*Specifically for the heat maps*)


ldplot[deldegcoords_,zcoords_]:=ldplot[deldegcoords,zcoords,{"n","#[delays]"},All,"TemperatureMap"];


(* ::Subsubsection:: *)
(*ClusterPlot of a 2D array*)


getclusterplot[zmat_, nclusters_] := Module[{zmatM, zmatN, zmatindices, zmatclusters},
   {zmatM, zmatN} = Dimensions@zmat;
   zmatindices = Table[{i, j}, {i, 1, zmatM}, {j, 1, zmatN}];
   zmatclusters = FindClusters[Flatten[zmat] -> Flatten[zmatindices, 1], nclusters];
   MatrixPlot[zmatindices /. Flatten@MapThread[(*Points in the same cluster are assigned the same value from 1 to nclusters *)Thread[#1 -> #2] &, {zmatclusters, Range[nclusters]}]]
   ];
 


(* ::Text:: *)
(*Something more suited to crunchcoords*)


getclusterplot[zmat_, zmatindices_, nclusters_] := Module[{zmatclusters, clusteredzmat},
   zmatclusters = FindClusters[Flatten[zmat] -> zmatindices, nclusters];
   clusteredzmat = zmatindices /. Flatten@MapThread[(*Points in the same cluster are assigned the same value from 1 to nclusters *)Thread[#1 -> #2] &, {zmatclusters, Range[nclusters]}];
   ldplot[zmatindices, clusteredzmat]
   ];


(* ::Subsubsection:: *)
(*Generate heat-maps of numerical checks from vals*)


vals2plotopsVScases[vals_,crunchcoords_,ncases_,fun_]:=Module[{casewiserates,opsVScases,funnyopsVScases,plotopsVScases},
casewiserates = (* Case, IC,del,deg*)Table[ Map[(#[ratequads])[[All,All,i]]&,vals],{i,ncases}];opsVScases = Table[casewiserates[[j,All,All,All,i]],{i,3}(* PDF, rperror, respercent *),{j,ncases}(* Vanilla, MS, MS + Pres, MS + Delay *)];
funnyopsVScases = Map[makecratesfunny[fun,#&,#]&,Rest[opsVScases],{2}];
plotopsVScases = Map[ldplot[crunchcoords,splog10@#]&,funnyopsVScases ,{2}];
MatrixForm[plotopsVScases\[Transpose]]
]


vals2plotopsVScases[vals_,crunchcoords_,ncases_]:= vals2plotopsVScases[vals,crunchcoords,ncases,Quantile[#,95/100]&]
