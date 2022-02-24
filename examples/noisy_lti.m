(* ::Package:: *)

(* ::Chapter:: *)
(*Initialization*)


(* ::Subsubsection:: *)
(*Run scripts containing functions*)


Get[FileNameJoin[{DirectoryName[NotebookFileName[]], "makefile.m"}]];



Get[FileNameJoin[{Nest[DirectoryName, NotebookFileName[],2],"src","function_forge.m"}] ];


(* ::Subsubsection::Closed:: *)
(*Setup Associations (structs in Common) to collect input parameters*)


{trajgrub,liftgrub,crunchgrub,plotgrub} = Table[<||>,{i,4}];


(* ::Subsubsection::Closed:: *)
(*Generate savefile name*)


AssociateTo[crunchgrub,savefile ->FileBaseName[NotebookFileName[]](* Name savefile after notebook *) ];
crunchgrub[plotdir] = FileNameJoin[{FileNameJoin[Drop[FileNameSplit[NotebookFileName[]],-2]],"plots"}];


(* ::Chapter:: *)
(*Computation*)


(* ::Section:: *)
(*Temporal parameters*)


(* ::Subsection:: *)
(*Define the scope of study*)


(* ::Item:: *)
(*Delays: mindelays \[LongRightArrow] maxdelays*)


(* ::Item:: *)
(*DMD model order (n) : minn \[LongRightArrow] maxn*)


AssociateTo[plotgrub,{testdelays -> Range[0,25](*Range[6,25]*)(*{6,25,50,100,200,400}*), testdegs-> Range[2,25]}];


(* ::Subsubsection::Closed:: *)
(*Modify some structs with the given*)


crunchgrub = Merge[{crunchgrub, plotgrub}, First];
crunchgrub[maxdelays] = Max@ crunchgrub[testdelays];



(* ::Subsection:: *)
(*DS specs*)


(* ::Item:: *)
(*LTI*)


(* ::Text:: *)
(*discevals: Discrete time eigenvalues of the LTI system *)


(* ::Subitem:: *)
(*LTI_{1}: Exp[I*2\[Pi]/7*Range[0,6]]*)


(* ::Subsubitem:: *)
(*LTI_{1.5} : Exp[I*2\[Pi]/7*Range[0,6]]*Exp[I*2\[Pi]/14] *)


(* ::Subitem:: *)
(*LTI_{2}: RandomReal[{0.8,1.2},7]*Exp[I*RandomReal[{-\[Pi],\[Pi]},7]]*)


(* ::Item:: *)
(*VDP : x1'=x2,x2'=\epsilon(1-x1^2)x2-x1 .*)


(*AssociateTo[trajgrub, discevals ->RandomReal[{0.8,1.2},7]*Exp[I*RandomReal[{-\[Pi],\[Pi]},7]]  (*Exp[I*(2\[Pi])/7*Range[0,6]]*)(**Exp[\[ImaginaryI]*(2\[Pi])/14] *)]; *)


(* ::Text:: *)
(*UnitDisc + UnitCircle*)


trajgrub[discevals] = Join[RandomReal[{1,1},3]*Exp[I*RandomReal[{-\[Pi],\[Pi]},3]],RandomReal[{1,1},4]*Exp[I*RandomReal[{-\[Pi],\[Pi]},4]]];


(* ::Subsubsection::Closed:: *)
(*ssdim:= Dimension of the LTI = "r"*)


ssdim = Length@trajgrub[discevals];


(* ::Subsection:: *)
(*Dictionary specs*)


(* ::Item:: *)
(*cmatseed*)


(* ::Item:: *)
(*crank: Rank of \tilde{C}*)


(* ::Item:: *)
(*crows: #(Rows) of  \tilde{C} = "m"*)


(* ::Item:: *)
(*nprojs : #(Non-zero components) for initial conditions*)


(* ::Subitem:: *)
(*nprojs = ssdim \[DoubleLongRightArrow] C and \tilde{C} have the same number of columns ("r")*)


(* ::Subitem:: *)
(*In addition,  crank = ssdim & crows >= ssdim ensures that "C" has full column rank*)


AssociateTo[liftgrub,{cmatseed-> {},crank ->(*ssdim *)1, crows -> 1 (**)  , nprojs -> ssdim (* \[LessEqual] ssdim ...tis simply the number of non-zero entires in the IC*)}];

crunchgrub[meff] = 2*ssdim;


(* ::Subsection:: *)
(*Sampling time*)


(* ::Item:: *)
(*tinit: Initial time for ODE*)


(* ::Subitem:: *)
(*Irrelevant ATM as all DS are autonomous*)


(* ::Item:: *)
(*tsamp: Time step to sample continuous flow*)


AssociateTo[trajgrub,{tinit -> 0, tsamp -> 1}];


(* ::Subsubsection::Closed:: *)
(*Generate integrating parameters*)


AssociateTo[trajgrub,{maxdelays ->  Max@crunchgrub[testdelays], maxn-> Max@crunchgrub[testdegs]}];


simsteps=2*trajgrub[maxdelays]+trajgrub[maxn]+1+crunchgrub[meff];


tfinn = trajgrub[tinit] + trajgrub[tsamp]*simsteps;
sampackage = {trajgrub[tinit],trajgrub[tsamp],simsteps-1};


(* ::Subsubsection::Closed:: *)
(*Find continuous time version of LTI*)


contevals = 1/trajgrub[tsamp]*Log[trajgrub[discevals]];
generatingAmat = eval2amat[contevals];


(* ::Subsection:: *)
(*Mean subtraction data*)


(* ::Item:: *)
(*rate2sub : "Number" whose mean is being removed (\[Lambda])*)


(* ::Subitem:: *)
(* 1 for all examples*)


(* ::Item:: *)
(*nfunda : Parameters used to decide when Case 1 can occur, and accordingly change the reference for \[Rho]_{Subset} and \[Delta]_{Trivial} using  getevalpartitions (See Table 6.1)*)


(* ::Subitem:: *)
(*This is of import purely when studying LTI systems whose spectrum is known*)


(* ::Subsubitem:: *)
(*Subscript[LTI, 1]: 6*)


(* ::Subsubitem:: *)
(*Subscript[LTI, 1.5]: 13*)


(* ::Subitem:: *)
(*All others : 2* Maximum value of "n" simulated*)


(* ::Item:: *)
(*truevals: True spectrum of the underlying system - Altered after all computations to check effectiveness of Algorithm 6.1*)


(* ::Subitem:: *)
(*LTI: trajgrub[discevals]*)


(* ::Subitem:: *)
(*Others: {}*)


AssociateTo[liftgrub,{rate2sub -> 1,nfunda-> (* Case 2 always *) 2*trajgrub[maxn] (* Impute the number you've rigged the system to have *)(*6*)(* 13*)}];
AssociateTo[liftgrub,truevals-> trajgrub[discevals](* {} if you don't know what it should be, in which case we don't know what mean subtraction does *)]; 


(* ::Subsubsection:: *)
(*Choice of DS*)


(* ::Item:: *)
(*vfield : Continuous time DS*)


(* ::Subitem:: *)
(*I/P: vector-field[time, coordinate vector ]*)


(* ::Item:: *)
(*chosenputty: Coordinate transform on ODE state into dictionary inputs*)


(* ::Subitem:: *)
(*Redundant for examples in draft*)


(* Create a separate one for each DS*)
locallti[t_(* Time *),r_(* Vector of coordinates*)]:=generatingAmat . r;
(* Assign vfield*)
AssociateTo[trajgrub,vfield ->  locallti];
(* chosenputty: Use to convert coordinates from the ODE version *)
AssociateTo[trajgrub, chosenputty -> (#&) ];


(* ::Subsection:: *)
(*Tolerances *)


(* ::Item:: *)
(*sigtols : reltol,abstol*)


(* ::Subitem:: *)
(*Former used in sigthresholder to implement the relative error criterion in SVD*)


(* ::Item:: *)
(*restols*)


(* ::Subitem:: *)
(*[[1]]:  Remove \[Lambda] from a chosen set in removelambfromset*)


(* ::Subitem:: *)
(*[[2]]:  Decide if a coefficient is 0 or not in a c-vector via maketrim*)


(* Tolerance *)
AssociateTo[crunchgrub,{sigtols ->  {10^-8,10^-12} (* Singular values *),
restols ->{10^-8,10^-12}(* greaterabsrelcheck *)}];
(**)


(* ::Subsubsection::Closed:: *)
(*Parameters for cte2wt *)


(* ::Item:: *)
(*coredmd computes a reduced SVD of "X" so that the condition number is below a threshold specified by restol*)


(* ::Subitem:: *)
(*( Condition number of the reduced "X" matrix )^cpow*)


(* ::Subitem:: *)
(*( Largest singular value dropped in reduction of X )^tepow*)


AssociateTo[crunchgrub,{cpow-> 1,tepow-> 1}];


(* ::Subsubsection:: *)
(*Noise parameters:*)


trajgrub[basicdist] = UniformDistribution[{-Sqrt[3],Sqrt[3]}] (*1*); (* Zero mean and SD = 1 *)
trajgrub[noiseSD]= 5;


trajgrub[covmat]= trajgrub[noiseSD]^2*IdentityMatrix[liftgrub[crows]];


(* ::PageBreak:: *)
(**)


(* ::Subsection:: *)
(*#[realizations] to analyse*)


(* ::Item:: *)
(*nICs = #[trajectories] to analyse*)


nICs =50;


(* ::Subsubsection::Closed:: *)
(*Generate a list of trajectories*)


listoICs= Table[getic[ssdim,liftgrub[nprojs]],{i,nICs}];
listotseries = Map[getimeseries[#,trajgrub[vfield],trajgrub[tinit],tfinn,sampackage,trajgrub[chosenputty]]&,listoICs];


listopnoise=Table[Transpose@getIIDnoise[{liftgrub[crows],simsteps+1},trajgrub[basicdist],trajgrub[covmat]],{i,nICs}];


(* ::Subsubsection::Closed:: *)
(*Initialize colour-scheme for plots*)


basicolourlist = Array[Hue[#]&,Length@crunchgrub[testdelays],{0,0.7 (* The end of the spectrum before it starts repeating *)}];


(* ::Section:: *)
(*Analyse each trajectory, for all possible hyper-parameters.*)


(* ::Input:: *)
(**)
(*trajgrub[noisyQ]= True;*)


(* ::Input:: *)
(*Print[AbsoluteTime[]];*)
(*vals=ParallelTable[*)
(*(trajgrub[rawics]=listotseries[[i]];*)
(*trajgrub[opnoise]=listopnoise[[i]];*)
(*meansuboneshot[trajgrub,liftgrub,crunchgrub,<||>]*)
(*(*Last bit changes when we have already done one pass,and wish to update our numerical checks*)),*)
(*{i,nICs}];*)
(*Print[AbsoluteTime[]];*)


(* ::Subsubsection:: *)
(*Save progress*)


(* ::Input:: *)
(*DumpSave[crunchgrub[savefile],{vals,crunchgrub,trajgrub,liftgrub,basicolourlist,listoICs,listotseries,plotgrub}];*)


(* ::Chapter:: *)
(*Restore-point #1*)


(* ::Input:: *)
(*(*Get[crunchgrub[savefile]];*)*)


(* ::PageBreak:: *)
(**)


(* ::Chapter:: *)
(*Post-processing*)


(* ::Subsubsection::Closed:: *)
(*Plot travails*)


(* ::Text:: *)
(*Legend for colour scheme*)


(* ::Input:: *)
(*delayscolored = BarLegend[{basicolourlist,Through[{Min,Max}[crunchgrub[testdelays]]]},crunchgrub[testdelays],LegendLabel->"#[delays]",LabelStyle->{Directive[Black,15]},LegendMarkerSize->{250}];*)


(* ::Text:: *)
(*Locate position of "r" within degrees("n") being tested*)


(* ::Input:: *)
(*findnprojsintestdegs = (Flatten[Position[crunchgrub[testdegs],liftgrub[nprojs]]])[[1]];*)


(* ::Input:: *)
(*(* Process you data before plotting - Makes labelling more intuitive *)*)
(*basicBWplot[ensembledat_,colourlist_,{vertmin_,vertmax_},vertcolour_]:=Module[{tplots,xlabelspacing = 3,sparseXlabels},*)
(*(* Plot only a subset of xlabels to avoid cluttering for larger fonts *)*)
(*(*Replace those whose index %xlabelspacing != 1 with None *)*)
(*sparseXlabels = sparsifyChartLabels[crunchgrub[testdegs],xlabelspacing];*)
(*tplots= MapThread[*)
(*BoxWhiskerChart[#1,{{"Whiskers",Transparent},{"Fences",Transparent},{"MedianMarker","o",#2}},ChartStyle->#2,ChartLabels-> sparseXlabels,ImageSize->Large,PlotRange->All,LabelStyle->Directive[Black, Medium],Joined->True,BarSpacing ->0.6 (* Allows the Join of Medians to be discerned, if reader is interested *),  PlotTheme->"Scientific"]&,*)
(*{ensembledat,colourlist}];*)
(*( *)
(*Show[##,ParametricPlot[{findnprojsintestdegs,t},{t,vertmin,vertmax},PlotStyle->Directive[vertcolour,Dashed](*Hue[0.8]*),PlotRange->All,ImageSize->Large,LabelStyle->Directive[Black,Medium]],PlotRange-> All]&*)
(*)@@(tplots)*)
(*];*)
(**)
(*stdBWplot[ensembledat_,colourlist_,colorlegend_,xlabel_,ylabel_,vertbounds_,vertcolour_]:=Module[{tplots,transitplot},transitplot=basicBWplot[ensembledat,colourlist,vertbounds,vertcolour];*)
(*(*Put the correct legend*)(Legended[#,colorlegend]&)@(*On this bunch of plots that have been overlaid*)Show[transitplot,FrameLabel->{{ylabel,None},{xlabel,None}},PlotLabel->None,LabelStyle->{Directive[Black,20]}]];*)


(* ::Input:: *)
(*kmdplots[key_, fun_, vals_, funname_] := Module[{funnyvals, plots},*)
(*   (* Cases, Del, Deg, ICs *)*)
(*   funnyvals = Transpose[*)
(*     Map[*)
(*      Map[fun, #[key], {3}] &,*)
(*      vals],*)
(*     4 <-> 1*)
(*     ];*)
(*   plots = *)
(*    Map[stdBWplot[#, basicolourlist, delayscolored, "n", *)
(*       funname, {0, 1}, Hue[0.8]] &, funnyvals];*)
(*   plots*)
(*   ];  *)


(* ::Chapter:: *)
(*Noisy plots*)


(* ::Input:: *)
(*plotgrub[casestrings] = {"noise_free","noisy_Companion","TLS","noise_resist_DMD"};*)


(* ::Subsection:: *)
(*Performance of varying delays for each case*)


(* ::Input:: *)
(*tplots = kmdplots[ratequads,kmdQuality,vals,"KMDQuality"];*)


(* ::Input:: *)
(*savetheseplots[tplots,nametheseplots[#,Map[StringJoin["KMDQuality_",#]&,plotgrub[casestrings]]]&,"png"];*)


(* ::Subsection:: *)
(*Best delay for all cases*)


(* ::Input:: *)
(*algocolours = {Green, Red, Cyan, Blue};*)
(*algoscoloured = LineLegend[algocolours, {"Companion", "Noisy Companion", "TLS", "Noise-resist Companion"}];*)


(* ::Input:: *)
(*Block[{keep = {1,2,3,4},delayindex =1},stdBWplot[kmdQuals[[keep,delayindex]],algocolours[[keep]],algoscoloured, "n","KMDQuality",{0,1},Hue[0.8]]]*)


(* ::Chapter::Closed:: *)
(*Restore-point #2*)


(* ::Text:: *)
(*Complete save-point*)


(* ::Input:: *)
(*(*Get[crunchgrub[savefile]];*)*)
