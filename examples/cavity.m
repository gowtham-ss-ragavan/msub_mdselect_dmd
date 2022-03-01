(* ::Package:: *)

(* ::Chapter:: *)
(*Initialization*)


(* ::Subsubsection::Closed:: *)
(*Run scripts containing functions*)


Get[FileNameJoin[{DirectoryName[NotebookFileName[]], "makefile.m"}]];



Get[FileNameJoin[{Nest[DirectoryName, NotebookFileName[],2],"src","function_forge.m"}] ];


(* ::Subsubsection::Closed:: *)
(*Setup Associations (structs in Common) to collect input parameters*)


{trajgrub,liftgrub,crunchgrub,plotgrub} = Table[<||>,{i,4}];


(* ::Subsection:: *)
(*Savefile // Choice of velocity field*)


(* ::Item:: *)
(*Valid options : 13k, 16k, 20k, 30k*)


AssociateTo[trajgrub, vfield -> "30k"]; 


(* ::Subsubsection::Closed:: *)
(*Generate savefile name*)


(* ::Input:: *)
(*AssociateTo[crunchgrub, savefile -> StringJoin[{FileBaseName[NotebookFileName[]] , "_", trajgrub[vfield]}]];*)
(**)
(*crunchgrub[plotdir] = FileNameJoin[{FileNameJoin[Drop[FileNameSplit[NotebookFileName[]],-2]],"plots"}];*)


(* ::Chapter:: *)
(*Computations*)


(* ::Section:: *)
(*Temporal parameters*)


(* ::Subsection::Closed:: *)
(*Define the scope of study*)


(* ::Item:: *)
(*Delays: mindelays \[LongRightArrow] maxdelays*)


(* ::Item:: *)
(*DMD model order (n) : minn \[LongRightArrow] maxn*)


AssociateTo[plotgrub,{testdelays -> Range[0,28](*{6,25,50,100,200,400}*), testdegs-> Range[2,25]}];


(* ::Subsubsection:: *)
(*Modify some structs with the given*)


crunchgrub = Merge[{crunchgrub, plotgrub}, First];
crunchgrub[maxdelays] = Max@ crunchgrub[testdelays];



(* ::Subsection::Closed:: *)
(*DS specs*)


(* ::Input:: *)
(*AssociateTo[trajgrub, chosenputty -> (#&) ];*)


(* ::Subsubsection::Closed:: *)
(*Get addresses and load the cavity data*)


(* ::Input:: *)
(*AssociateTo[trajgrub, datadir->FileNameJoin[{ParentDirectory[DirectoryName[NotebookFileName[]]],"data"}]]; *)


(* ::Input:: *)
(*AssociateTo[trajgrub,prunedata-> (FileNames[__~~trajgrub[vfield]~~__,FileNameJoin[{trajgrub[datadir],"pruned"}]])[[1]]];*)


(* ::Input:: *)
(*Get[trajgrub[prunedata]];*)
(**)


(* ::Input:: *)
(*plotgrub[casestrings] = {"vanilla","ms","mspres","msdel"};*)


(* ::Subsubsection::Closed:: *)
(*ssdim:= #[grid-points] that sample*)


(* ::Input:: *)
(*{ssdim,rawsnapcount}=Dimensions[velocityfield[cavityPsi]];*)


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Sampling time*)


(* ::Item:: *)
(*tinit: Initial time for ODE*)


(* ::Subitem:: *)
(*Irrelevant ATM as all DS are autonomous*)


(* ::Item:: *)
(*tsamp: Time step to sample continuous flow*)


(* ::Input:: *)
(*AssociateTo[trajgrub,{tinit -> 0, tsamp -> 0.05}];*)


(* ::Subsubsection::Closed:: *)
(*Generate integrating parameters*)


AssociateTo[trajgrub,{maxdelays ->  Max@crunchgrub[testdelays], maxn-> Max@crunchgrub[testdegs]}];


simsteps=2*trajgrub[maxdelays]+trajgrub[maxn]+1+crunchgrub[meff];


(* ::Subsection::Closed:: *)
(*Mean subtraction data*)


(* ::Item:: *)
(*rate2sub : "Number" whose mean is being removed (\[Lambda])*)


(* ::Subitem:: *)
(* 1 for all examples*)


(* ::Item:: *)
(*nfunda : Parameters used to decide when Case 1 can occur, and accordingly change the reference for  \[Rho]_{Subset} and \[Delta]_{Trivial} using  getevalpartitions (See Table 6.1)*)


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


(* ::Input:: *)
(*AssociateTo[liftgrub,{rate2sub -> 1,nfunda-> (* Case 2 always *) 2*trajgrub[maxn] (* Impute the number you've rigged the system to have *)(*6*)(* 13*)}];*)


(* ::Subsection::Closed:: *)
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


(* ::Input:: *)
(*(* Tolerance *)*)
(*AssociateTo[crunchgrub,{sigtols ->  {10^-8,10^-12} (* Singular values *),*)
(*restols ->{10^-8,10^-12}(* greaterabsrelcheck *)}];*)
(*(**)*)


(* ::Subsubsection::Closed:: *)
(*Parameters for cte2wt *)


(* ::Item:: *)
(*coredmd computes a reduced SVD of "X" so that the condition number is below a threshold specified by restol*)


(* ::Subitem:: *)
(*( Condition number of the reduced "X" matrix )^cpow*)


(* ::Subitem:: *)
(*( Largest singular value dropped in reduction of X )^tepow*)


(* ::Input:: *)
(*AssociateTo[crunchgrub,{cpow-> 1,tepow-> 1}];*)


(* ::Subsubsection::Closed:: *)
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


(* ::Input:: *)
(*nICs = 50;*)


(* ::Subsubsection:: *)
(*Generate a list of trajectories (special to cavity)*)


(* ::Input:: *)
(*choices4IC =Range[Floor[ 0.8*rawsnapcount]] (* Range[1000] *);*)
(*(* If choices4IC is too small, will abort*)*)
(*If[(Max[choices4IC] + simsteps) > rawsnapcount,Print["Insufficient data : make choices4IC OR simsteps smaller"]; Quit[]];*)
(*(* If nICs is too greedy, will Quit too *)*)
(*If[ nICs > Length[choices4IC], Print[" Insufficient data: Reduce nICs OR Increase choices4IC"];Quit[]];*)
(**)


(* ::Input:: *)
(*(* So, we only store like the start points *)*)
(*listoICs= RandomSample[choices4IC,nICs];*)
(*listotseries = Map[((velocityfield[cavityPsi])[[All,Range[0,simsteps]+#]])\[Transpose](* Transpose coz the cavityPsi is in snapshot form *)&,listoICs];*)


(* ::Input:: *)
(*listopnoise=Table[Transpose@getIIDnoise[{liftgrub[crows],simsteps+1},trajgrub[basicdist],trajgrub[covmat]],{i,nICs}];*)


(* ::Subsubsection:: *)
(*Initialize colour-scheme for plots*)


(* ::Input:: *)
(*basicolourlist = Array[Hue[#]&,trajgrub[maxdelays]+1,{0,0.7 (* The end of the spectrum before it starts repeating *)}];*)


(* ::Input:: *)
(*(*Abs[trajgrub[discevals]]*)*)


(* ::Section:: *)
(*Analyse each trajectory, for all possible hyper-parameters.*)


(* ::Input:: *)
(*liftgrub[truevals] = {};*)
(*trajgrub[noisyQ]= False;*)


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
(*(*Get[crunchgrub[savefile]];*)
(*(* Load the velocityfield data *)*)
(*Get[trajgrub[prunedata]];*)
(*(* Generate your trajectories *)*)
(*nICs = Length[listoICs];*)
(*listotseries = Map[((velocityfield[cavityPsi])[[All,Range[0,simsteps]+#]])\[Transpose](* Trasnpose coz the cavityPsi is in snapshot form *)&,listoICs];*)
(*listopnoise=Table[Transpose@getIIDnoise[{liftgrub[crows],simsteps+1},trajgrub[basicdist],trajgrub[covmat]],{i,nICs}];*)
(*(* UNdo all changes *)*)
(*If[ValueQ[oldvals],liftgrub[truevals]={};*)
(*KeyDropFrom[crunchgrub,{rmin,delaymin,chosendeg}];vals = oldvals];*)*)


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


(* ::Item:: *)
(*No use for VDP and cavity - > IS this the best inert value it can be given ?*)


(* ::Input:: *)
(*findnprojsintestdegs = (crunchgrub[testdegs])[[1]] ;*)


(* ::Section:: *)
(*DMD-DFT transition*)


(* ::Subsubsection:: *)
(*Parse data into a matrix form*)


(* ::Text:: *)
(*Compute \delta_{rel}^{\mu}*)


(* ::Input:: *)
(*ensembledat = Transpose[Map[ Map[obs4c,(#[cquads])[[All,All,2]],{2}]&,vals],{3,1,2}];*)


(* ::Text:: *)
(*Collect \[Kappa](X),\sigma_{tail}(X)*)


(* ::Input:: *)
(*ensembleconditiondat = Transpose[Map[ ((#[equads])[[All,All,2,1]])&,vals],{3,1,2}];*)
(*ensemblechoppeddat = Transpose[Map[ ((#[equads])[[All,All,2,3 (* {cond,condred,tail}*)]])&,vals],{3,1,2}];*)


(* ::Subsection:: *)
(* Deviation of mean-subtracted DMD from DFT Vs n*)


(* ::Item:: *)
(*Find \hat{r}_{min}, \hat{r}_{max} so that when at least n_{max}-1 delays are taken,  the jump in the indicator  occurs for  \hat{r}_{min} <=  n <= \hat{r}_{max}*)


(* ::Input:: *)
(*dmddftdeviationplot = stdBWplot[crunchgrub[testdegs],splog10@ensembledat,basicolourlist,delayscolored,HoldForm[n],HoldForm[Subscript[Log, 10][\!\(\*SubsuperscriptBox[\(\[Delta]\), \(rel\), \(\[Mu]\)]\)]],{-17,1},Transparent]*)


(* ::Input:: *)
(*savetheseplots[{dmddftdeviationplot },nametheseplots[#,{"dmddftdeviation"}]&,"png"];*)


(* ::Subsubsection:: *)
(*Tail of SVD lost in truncation Vs n*)


(* ::Input:: *)
(*stdBWplot[crunchgrub[testdegs],splog10@ensemblechoppeddat,basicolourlist,delayscolored,HoldForm[n],HoldForm[Subscript[Log, 10][Subscript[\[Sigma], Tail]]],{-17,1},Transparent]*)


(* ::Section:: *)
(*Eigenvalues from multiple time traces : Algorithm 6.1*)


(* ::Subsection:: *)
(*Determining rmin*)


(* ::Subsubsection:: *)
(*Log10[How well is the y-th estimate contained in the x-th estimate] *)


(* ::Input:: *)
(*howtochooseyourrplot =getetvdepwrtcdeg[vals,crunchgrub]*)


savetheseplots[{howtochooseyourrplot},nametheseplots[#,{"rpyramid"}]&,"png"];


(* ::Subsection:: *)
(*Choose  rmin from the earlier diagnostics*)


(* ::Item:: *)
(*Highest y-coordinate, within \hat{r}_{min}: \hat{r}_{max}, with a row of low values*)


(* ::Input:: *)
(*AssociateTo[crunchgrub,{rmin-> Input["Please enter your guess of rmin ('hat{r}')"](* 7*),chosendeg->  Input["Please enter hat{n} to sample hat{r} from "]}];*)


(* ::Subsubsection:: *)
(*Pick data with sufficient delays and extract their common eigenvalues*)


(* ::Input:: *)
(*AssociateTo[crunchgrub,{delaymin-> crunchgrub[chosendeg]-1}];*)
(*valswithgoodelays=keepthemgoodelays[onlythequads[vals],crunchgrub];*)


(* ::Input:: *)
(*{prunedrootqfs,estruevals}=sneakyestimatetruevals[valswithgoodelays,crunchgrub,"verbose"];*)
(*Clear[valswithgoodelays];*)


(* ::Subsubsection:: *)
(*Minutia*)


(* ::Input:: *)
(*oldvals = vals;*)
(*ncases=4;*)


(* ::Subsubsection:: *)
(*Generate coordinates for heatmaps*)


(* ::Input:: *)
(*plotcoords = Flatten[Outer[{##}&,plotgrub[testdelays],plotgrub[testdegs]],1];*)
(*crunchcoords = Flatten[Outer[{##}&,crunchgrub[testdelays],crunchgrub[testdegs]],1];*)


(* ::Section:: *)
(*Theorem checks*)


(* ::Item:: *)
(*Each case generates a 4 x 2 table of heat-maps, each plot describing the variation of an averaged quantity with respect to #[delays] and n*)


(* ::Subitem:: *)
(*Columns : \[Rho]_{Subset}, \[Delta]_{Trivial}*)


(* ::Subitem:: *)
(*Rows: Theorems for Vanilla, Mean-subtracted, mspres and msdel*)


(* ::Item:: *)
(*All theorems predict that beyond a critical value of n determined by r (the dimension of the underlying Koopman invariant subspace), all heat maps* should exhibit low values*)


(* ::Subitem:: *)
(*The second row is a little subtle and better explained in the manuscript*)


(* ::Subsection:: *)
(*Estimated eigenvalues: \!\(\*OverscriptBox[\(\[Nu]\), \(^\)]\)*)


(* ::Input:: *)
(*liftgrub[truevals] = estruevals;*)
(*(* We want the vertical line at the spot where we hit the number of eigenvalues used as the truth in the current analysis *)*)
(*findnprojsintestdegs = (Flatten[Position[crunchgrub[testdegs],Length[liftgrub[truevals]]]])[[1]];*)
(**)
(*vals = ParallelMap[*)
(*ms1shot4trajvariations[trajgrub, liftgrub, crunchgrub, #] &, *)
(*{listotseries, listopnoise, oldvals}\[Transpose]*)
(*];*)


(* ::Input:: *)
(*estimateplots = kmdplots[crunchgrub[testdegs],ratequads,kmdQuality,vals,"KMD-Quality",basicolourlist,delayscolored]*)


(* ::Input:: *)
(*savetheseplots[estimateplots,nametheseplots[#,prefixstrlist["estimatevals_",plotgrub[casestrings]]]&,"png"];*)


(* ::Subsubsection:: *)
(*Save again*)


(* ::Input:: *)
(*DumpSave[crunchgrub[savefile],{oldvals,vals,crunchgrub,trajgrub,liftgrub,basicolourlist,listoICs,plotgrub,simsteps,crunchcoords,ncases}];*)


(* ::PageBreak:: *)
(**)


(* ::Chapter:: *)
(*Restore-point #2*)


(* ::Text:: *)
(*Complete save-point*)


(* ::Input:: *)
(*(*Get[crunchgrub[savefile]];*)
(*(* Load the velocityfield data *)*)
(*Get[trajgrub[prunedata]];*)
(*(* Generate your trajectories *)*)
(*listotseries = Map[((velocityfield[cavityPsi])[[All,Range[0,simsteps]+#]])\[Transpose](* Trasnpose coz the cavityPsi is in snapshot form *)&,listoICs];*)
(*listopnoise=Table[Transpose@getIIDnoise[{liftgrub[crows],simsteps+1},trajgrub[basicdist],trajgrub[covmat]],{i,nICs}];*)
(**)*)


Sort[Abs[estruevals]]


(* ::Text:: *)
(*20k: {0.999255, 1.06807}*)


(* ::Text:: *)
(*16k:  {0.994866, 0.997453}*)


(* ::Text:: *)
(*13k: {0.965107, 0.965107, 0.991791, 0.991791, 0.999884}*)
