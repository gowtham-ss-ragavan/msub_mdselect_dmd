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


AssociateTo[trajgrub, vfield -> "13k"]; 


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


(* ::Subsection:: *)
(*Define the scope of study*)


(* ::Item:: *)
(*Delays: mindelays \[LongRightArrow] maxdelays*)


(* ::Item:: *)
(*DMD model order (n) : minn \[LongRightArrow] maxn*)


AssociateTo[plotgrub,{testdelays -> {0, 12, 25}(*Range[0,28]*)(*{6,25,50,100,200,400}*), testdegs-> Range[2,26]}];


(* ::Subsubsection:: *)
(*Modify some structs with the given*)


crunchgrub = Merge[{crunchgrub, plotgrub}, First];
crunchgrub[maxdelays] = Max@ crunchgrub[testdelays];



(* ::Subsection:: *)
(*DS specs*)


(* ::Input:: *)
(*AssociateTo[trajgrub, chosenputty -> (#&) ];*)


(* ::Subsubsection:: *)
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


(* ::Subsection:: *)
(*Sampling time*)


(* ::Item:: *)
(*tinit: Initial time for ODE*)


(* ::Subitem:: *)
(*Irrelevant ATM as all DS are autonomous*)


(* ::Item:: *)
(*tsamp: Time step to sample continuous flow*)


(* ::Input:: *)
(*AssociateTo[trajgrub,{tinit -> 0, tsamp -> 0.05}];*)


(* ::Subsubsection:: *)
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


basicolourlist = Array[Hue[#]&,Length@crunchgrub[testdelays],{0,0.7 (* The end of the spectrum before it starts repeating *)}];


(* ::Section::Closed:: *)
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
(*(*Get[crunchgrub[savefile]];*)*)


(* ::Chapter:: *)
(*Post-processing*)


(* ::Subsubsection:: *)
(*Plot travails*)


(* ::Text:: *)
(*Legend for colour scheme*)


(* ::Input:: *)
(*delayscolored = SwatchLegend[basicolourlist,crunchgrub[testdelays],LegendLabel->"Delays (d)",LabelStyle->{Directive[Black,15]}, LegendLayout->"ReversedColumn", LegendFunction->"Frame"];*)


(* ::Text:: *)
(*Locate position of "r" within degrees("n") being tested*)


(* ::Item:: *)
(*No use for VDP and cavity - > IS this the best inert value it can be given ?*)


(* ::Input:: *)
(*findnprojsintestdegs = (crunchgrub[testdegs])[[1]] ;*)


(* ::Section:: *)
(*DMD-DFT transition*)


(* ::Subsubsection::Closed:: *)
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
(*dmddftdeviationplot = stdBWplot[crunchgrub[testdegs],splog10@ensembledat,basicolourlist,delayscolored,"Model-order(\[Theta])","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[ Relative distance to DFT ]",{-17,1},Transparent]*)


(* ::Input:: *)
(*savetheseplots[{dmddftdeviationplot },nametheseplots[#,{"dmddftdeviation"}]&,"png"];*)


(* ::Subsubsection:: *)
(*Tail of SVD lost in truncation Vs n*)


(* ::Input:: *)
(*stdBWplot[crunchgrub[testdegs],splog10@ensemblechoppeddat,basicolourlist,delayscolored,HoldForm[n],HoldForm[Subscript[Log, 10][Subscript[\[Sigma], Tail]]],{-17,1},Transparent]*)
