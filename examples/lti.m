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


(* ::Subsubsection:: *)
(*Set LTI case*)


crunchgrub[lticase] = "3";


(* ::Subsubsection::Closed:: *)
(*Generate savefile name*)


(* ::Input:: *)
(*AssociateTo[crunchgrub,savefile ->StringJoin[FileBaseName[NotebookFileName[]],"_",crunchgrub[lticase]](* Name savefile after notebook + Case number of LTI system  *) ];*)
(*crunchgrub[plotdir] = FileNameJoin[{FileNameJoin[Drop[FileNameSplit[NotebookFileName[]],-2]],"plots"}];*)


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


AssociateTo[plotgrub,{testdelays -> {0, 6, 24}(*Join[{0,3,5,6},{13,20,27}]*)(*Range[0,28]*), testdegs-> Range[2,26]}];


(* ::Subsubsection::Closed:: *)
(*Modify some structs with the given*)


crunchgrub = Merge[{crunchgrub, plotgrub}, First];
crunchgrub[maxdelays] = Max@ crunchgrub[testdelays];



(* ::Subsection::Closed:: *)
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


(* ::Input:: *)
(*plotgrub[casestrings] = {"vanilla","ms","truesansrate","truewithrate"};*)


(* ::Input:: *)
(*trajgrub[discevals] = Switch[crunchgrub[lticase],*)
(*"1a", Exp[I*(2\[Pi])/7*Range[0,6]],*)
(*"1b",RandomReal[{1,1},7]*Exp[I*RandomReal[{-\[Pi],\[Pi]},7]],*)
(*"2",RandomReal[{0.8,1.2},7]*Exp[I*RandomReal[{-\[Pi],\[Pi]},7]],*)
(*"3",Join[RandomReal[{0.8,1},3]*Exp[I*RandomReal[{-\[Pi],\[Pi]},3]],RandomReal[{1,1},4]*Exp[I*RandomReal[{-\[Pi],\[Pi]},4]]]*)
(*];*)


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


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
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


AssociateTo[liftgrub,{rate2sub -> 1,nfunda-> 2*trajgrub[maxn](* Irrelevant at this point: Just a dummy variable until there is time for a cleanup*)}];
AssociateTo[liftgrub,truevals-> trajgrub[discevals](* {} if you don't know what it should be, in which case we don't know what mean subtraction does *)]; 


(* ::Subsubsection::Closed:: *)
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


nICs = 50;


(* ::Subsubsection::Closed:: *)
(*Generate a list of trajectories*)


listoICs= Table[getic[ssdim,liftgrub[nprojs]],{i,nICs}];
listotseries = Map[getimeseries[#,trajgrub[vfield],trajgrub[tinit],tfinn,sampackage,trajgrub[chosenputty]]&,listoICs];


listopnoise=Table[Transpose@getIIDnoise[{liftgrub[crows],simsteps+1},trajgrub[basicdist],trajgrub[covmat]],{i,nICs}];


(* ::Subsubsection:: *)
(*Initialize colour-scheme for plots*)


basicolourlist = Array[Hue[#]&,Length@crunchgrub[testdelays],{0,0.7 (* The end of the spectrum before it starts repeating *)}];


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
(*(*Get[crunchgrub[savefile]];*)*)


(* ::PageBreak:: *)
(**)


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


(* ::Input:: *)
(*findnprojsintestdegs = (Flatten[Position[crunchgrub[testdegs],liftgrub[nprojs]]])[[1]];*)


(* ::Chapter:: *)
(*Noise-free plots*)


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
(*dmddftdeviationplot = stdBWplot[crunchgrub[testdegs],splog10@ensembledat,basicolourlist,delayscolored,"Model-order(\[Theta])","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[ Relative distance to DFT ]",{-17,1},Hue[0.8]]*)


(* ::Text:: *)
(**)


(* ::Input:: *)
(*savetheseplots[{dmddftdeviationplot },nametheseplots[#,{"dmddftdeviation"}]&,"png"];*)


(* ::Subsubsection:: *)
(*Tail of SVD lost in truncation Vs n*)


(* ::Input:: *)
(*stdBWplot[crunchgrub[testdegs],splog10@ensemblechoppeddat,basicolourlist,delayscolored,HoldForm[n],HoldForm[Subscript[Log, 10][Subscript[\[Sigma], Tail]]],{-17,1},Hue[0.8]]*)


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
(*True eigenvalues: \sigma(\[CapitalLambda])*)


(* ::Input:: *)
(*liftgrub[truevals] = N@trajgrub[discevals];*)
(*(* Updating the vertical line for a different set of target eigenvalues *)*)
(*findnprojsintestdegs = (Flatten[Position[crunchgrub[testdegs],Length[liftgrub[truevals]]]])[[1]];*)
(**)
(*vals =ParallelMap[*)
(*ms1shot4trajvariations[trajgrub, liftgrub, crunchgrub, #] &, *)
(*{listotseries, listopnoise, vals}\[Transpose]*)
(*];*)


(* ::Input:: *)
(*truthplots = kmdplots[crunchgrub[testdegs],ratequads,kmdQuality,vals,"Model-order(\[Theta])","KMD-Quality",basicolourlist,delayscolored]*)


(* ::Input:: *)
(*(*savetheseplots[truthplots,nametheseplots[#,prefixstrlist["truevals_",plotgrub[casestrings]]]&,"png"];*)*)


(* ::Subsubsection:: *)
(*Save again*)


(* ::Input:: *)
(*DumpSave[crunchgrub[savefile],{vals,crunchgrub,trajgrub,liftgrub,basicolourlist,listoICs,plotgrub,simsteps,crunchcoords}];*)


(* ::PageBreak:: *)
(**)


(* ::Chapter:: *)
(*Restore-point #2*)


(* ::Text:: *)
(*Complete save-point*)


(* ::Input:: *)
(*(*Get[crunchgrub[savefile]];*)*)
