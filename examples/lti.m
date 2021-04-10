(* ::Package:: *)

(* ::Chapter:: *)
(*Initialization*)


(* ::Subsubsection::Closed:: *)
(*Run scripts containing functions*)


(* ::Input:: *)
(*Get[FileNameJoin[{DirectoryName[NotebookFileName[]],"makefile.m"}]];*)


(* ::Subsubsection::Closed:: *)
(*Setup Associations (structs in Common) to collect input parameters*)


(* ::Input:: *)
(*{trajgrub,liftgrub,crunchgrub,priorsgrub,testpriorsgrub,plotgrub} = Table[<||>,{i,6}];*)


(* ::Subsubsection::Closed:: *)
(*Generate savefile name*)


(* ::Input:: *)
(*AssociateTo[crunchgrub,savefile ->FileBaseName[NotebookFileName[]](* Name savefile after notebook *) ];*)


(* ::Chapter:: *)
(*Computation*)


(* ::Subchapter::Closed:: *)
(*First run *)


(* ::Text:: *)
(*No pre-existing save files*)


(* ::Section:: *)
(*Temporal parameters*)


(* ::Subsection:: *)
(*Define the scope of study*)


(* ::Item:: *)
(*Delays: mindelays \[LongRightArrow] maxdelays*)


(* ::Item:: *)
(*DMD model order (n) : minn \[LongRightArrow] maxn*)


(* ::Input:: *)
(*AssociateTo[plotgrub,{mindelays-> 0, minn-> 2,maxdelays-> 20,maxn-> 20}];*)


(* ::Subsubsection::Closed:: *)
(*Modify some structs with the given*)


(* ::Input:: *)
(*crunchgrub = Merge[{crunchgrub,plotgrub},First];*)
(**)
(*{crunchgrub[maxdelays],crunchgrub[maxn]} += {4, 4} (* Arbitrary offset *);*)
(**)
(*{plotgrub,crunchgrub}=Map[Merge[{#,<|testdelays-> Range[#[mindelays],#[maxdelays]],testdegs-> Range[#[minn],#[maxn]]|>},#[[1]]&]&,{plotgrub,crunchgrub}];*)


(* ::Subsection:: *)
(*DS specs*)


(* ::Item:: *)
(*LTI*)


(* ::Subitem:: *)
(*discevals: Discrete time eigenvalues of the LTI system *)


(* ::Subitem:: *)
(*Subscript[LTI, 1]: Exp[I*(2\[Pi])/7*Range[0,6]]*)


(* ::Subsubitem:: *)
(*Subscript[LTI, 1.5] : Exp[I*(2\[Pi])/7*Range[0,6]]*Exp[I*(2\[Pi])/14] *)


(* ::Subitem:: *)
(*Subscript[LTI, 2]: RandomReal[{0.8,1.2},7]*Exp[I*RandomReal[{-\[Pi],\[Pi]},7]]*)


(* ::Item:: *)
(*VDP : Subscript[\!\(\*OverscriptBox[\(x\), \(.\)]\), 1]=Subscript[x, 2], Subscript[\!\(\*OverscriptBox[\(x\), \(.\)]\), 2] = \[Epsilon](1-Subscript[x, 1]^2)Subscript[x, 2]-Subscript[x, 1]*)


(* ::Input:: *)
(**)
(*AssociateTo[trajgrub, discevals ->RandomReal[{0.8,1.2},7]*Exp[I*RandomReal[{-\[Pi],\[Pi]},7]]  (*Exp[\[ImaginaryI]*(2\[Pi])/7*Range[0,6]]*Exp[\[ImaginaryI]*(2\[Pi])/14] *)]; *)
(**)


(* ::Subsubsection::Closed:: *)
(*ssdim:= Dimension of the LTI = "r"*)


(* ::Input:: *)
(*ssdim = Length@trajgrub[discevals];*)


(* ::Subsection:: *)
(*Dictionary specs*)


(* ::Item:: *)
(*cmatseed*)


(* ::Item:: *)
(*crank: Rank of \!\(\*OverscriptBox[\(C\), \(~\)]\)*)


(* ::Item:: *)
(*crows: #(Rows) of \!\(\*OverscriptBox[\(C\), \(~\)]\) = "m"*)


(* ::Item:: *)
(*nprojs : #(Non-zero components) for initial conditions*)


(* ::Subitem:: *)
(*nprojs = ssdim \[DoubleLongRightArrow] C and \!\(\*OverscriptBox[\(C\), *)
(*OverscriptBox[\(~\), \(\ \)]]\) have the same number of columns ("r")*)


(* ::Subitem:: *)
(*In addition,  crank = ssdim & crows >= ssdim ensures that "C" has full column rank*)


(* ::Input:: *)
(*AssociateTo[liftgrub,{cmatseed-> {},crank ->ssdim (*1*), crows -> 2*ssdim (**)  , nprojs -> ssdim (* \[LessEqual] ssdim ...tis simply the number of non-zero entires in the IC*)}];*)
(**)


(* ::Subsection:: *)
(*Sampling time*)


(* ::Item:: *)
(*tinit: Initial time for ODE*)


(* ::Subitem:: *)
(*Irrelevant ATM as all DS are autonomous*)


(* ::Item:: *)
(*tsamp: Time step to sample continuous flow*)


(* ::Input:: *)
(*AssociateTo[trajgrub,{tinit -> 0, tsamp -> 1}];*)


(* ::Subsubsection::Closed:: *)
(*Generate integrating parameters*)


(* ::Input:: *)
(*AssociateTo[trajgrub,{maxdelays ->  Max@crunchgrub[testdelays], maxn-> Max@crunchgrub[testdegs]}];*)


(* ::Input:: *)
(**)
(*simsteps = trajgrub[maxdelays] + trajgrub[maxn] (*+ 1 (* Delay after mean sub*)*) + 1 (* liftedX *);*)
(*tfinn = trajgrub[tinit] + trajgrub[tsamp]*simsteps;*)
(**)


(* ::Input:: *)
(*sampackage = {trajgrub[tinit],trajgrub[tsamp],simsteps-1};*)


(* ::Subsubsection::Closed:: *)
(*Find continuous time version of LTI*)


(* ::Input:: *)
(*contevals = 1/trajgrub[tsamp]*Log[trajgrub[discevals]];*)
(*generatingAmat = eval2amat[contevals];*)


(* ::Subsection:: *)
(*Mean subtraction data*)


(* ::Item:: *)
(*rate2sub : "Number" whose mean is being removed (\[Lambda])*)


(* ::Subitem:: *)
(* 1 for all examples*)


(* ::Item:: *)
(*nfunda : Parameters used to decide when Case 1 can occur, and accordingly change the reference for Subscript[\[Rho], Subset] and Subscript[\[Delta], Trivial] using  getevalpartitions (See Table 6.1)*)


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
(*AssociateTo[liftgrub,truevals-> trajgrub[discevals](* {} if you don't know what it should be, in which case we don't know what mean subtraction does *)]; *)


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


(* ::Input:: *)
(*(* Create a separate one for each DS*)*)
(*locallti[t_(* Time *),r_(* Vector of coordinates*)]:=generatingAmat.r;*)
(*(* Assign vfield*)*)
(*AssociateTo[trajgrub,vfield ->  locallti];*)
(*(* chosenputty: Use to convert coordinates from the ODE version *)*)
(*AssociateTo[trajgrub, chosenputty -> (#&) ];*)


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


(* ::Input:: *)
(*(* Tolerance *)*)
(*AssociateTo[crunchgrub,{sigtols ->  {10^-8,10^-12} (* Singular values *),*)
(*restols ->{10^-8,10^-12}(* greaterabsrelcheck *)}];*)
(*(**)*)


(* ::Subsubsection:: *)
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
(*Flavours  (Redundant)*)




(* ::Input:: *)
(*(* AppendTo[userlist,(*Desired fheads*)];*)*)
(**)
(*AssociateTo[priorsgrub,userparams ->  {}];*)




(* ::Item:: *)
(*Should be a LIST of matrices !!!*)


(* ::Input:: *)
(*AssociateTo[priorsgrub,{(*(*  *)*)
(*invsamplesets \[Rule] {},*)
(*(*  *)*)
(*ipfuns(**)\[Rule] {},*)
(*opfuns(**)\[Rule] {},*)*)
(*rawuserlist ->  {}}];*)


(* ::Input:: *)
(*(*AssociateTo[testpriorsgrub,{*)
(*testinvsamplesets \[Rule] {},*)
(*testipfuns \[Rule] {},*)
(*testopfuns \[Rule] {}}];*)*)


(* ::Input:: *)
(*AssociateTo[liftgrub,{flavour ->  "polytri", flavourparams -> 1}];*)


(* ::PageBreak:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Box parametrization (Redundant)*)


(* ::Item:: *)
(*Utility in future versions*)


(* ::Item:: *)
(*Currently redundant due to *)


(* ::Subitem:: *)
(*Choice of flavour as "polytri"*)


(* ::Subitem:: *)
(*Definition of basislift deployed *)


(* ::Input:: *)
(*boxranges = {{-1,1},{-1,1}};*)
(*centercounts = {{35,35},{9,9}(* Must be the same as the first to ensure that the scaling is correctly done for trigtri *)};*)
(*nsampsperbox = {5,5};*)
(*AssociateTo[liftgrub,boxflavourparams ->  {boxranges,centercounts[[-1]]}];*)


(* ::Subsection:: *)
(*#[realizations] to analyse*)


(* ::Item:: *)
(*nICs = #[trajectories] to analyse*)


(* ::Input:: *)
(*nICs =15;*)


(* ::Subsubsection::Closed:: *)
(*Generate a list of trajectories*)


(* ::Input:: *)
(*listoICs= Table[getic[ssdim,liftgrub[nprojs]],{i,nICs}];*)
(*listotseries = Map[getimeseries[#,trajgrub[vfield],trajgrub[tinit],tfinn,sampackage,trajgrub[chosenputty]]&,listoICs];*)


(* ::Subsubsection::Closed:: *)
(*Initialize colour-scheme for plots*)


(* ::Input:: *)
(*basicolourlist = Array[Hue[#]&,trajgrub[maxdelays]+1,{0,0.7 (* The end of the spectrum before it starts repeating *)}];*)


(* ::Input:: *)
(*(*Abs[trajgrub[discevals]]*)*)


(* ::Section:: *)
(*First computing run*)


(* ::Subsection::Closed:: *)
(*Initial analysis without truevals*)


(* ::Input:: *)
(*liftgrub[truevals] = {};*)


(* ::Input:: *)
(*Print[AbsoluteTime[]];*)
(*vals =ParallelMap[(trajgrub[rawics] = #;meansuboneshot[trajgrub,liftgrub,priorsgrub,testpriorsgrub,crunchgrub,<||>(* To show that we do not have any prior knowledge of the system: This holds regardless of whether we know liftgrub[truevals] *)])&,listotseries];*)
(*Print[AbsoluteTime[]];*)


(* ::Input:: *)
(*DumpSave[crunchgrub[savefile],{vals,crunchgrub,trajgrub,liftgrub,priorsgrub,testpriorsgrub,basicolourlist,listoICs,listotseries,plotgrub}];*)


(* ::Subchapter:: *)
(*Later runs (i.e a savefile exists)*)


(* ::Input:: *)
(*(*Get[crunchgrub[savefile]];*)*)


(* ::PageBreak:: *)
(**)


(* ::Chapter:: *)
(*Post-processing*)


(* ::Section:: *)
(*Subscript[c, ms] transition*)


(* ::Subsubsection::Closed:: *)
(*Parse data into a matrix form*)


(* ::Text:: *)
(*Compute \!\(\*SubsuperscriptBox[\(\[Delta]\), \(rel\), \(\[Mu]\)]\)*)


(* ::Input:: *)
(*ensembledat = Transpose[Map[ Map[obs4c,(#[cquads])[[All,All,2]],{2}]&,vals],{3,1,2}];*)


(* ::Text:: *)
(*Collect \[Kappa](X),Subscript[\[Sigma], tail](X)*)


(* ::Input:: *)
(*ensembleconditiondat = Transpose[Map[ ((#[equads])[[All,All,2,1]])&,vals],{3,1,2}];*)
(*ensemblechoppeddat = Transpose[Map[ ((#[equads])[[All,All,2,3 (* {cond,condred,tail}*)]])&,vals],{3,1,2}];*)


(* ::Subsubsection::Closed:: *)
(*Plot travails*)


(* ::Text:: *)
(*Legend for colour scheme*)


(* ::Input:: *)
(*BarLegend[{basicolourlist,Through[{Min,Max}[crunchgrub[testdelays]]]},crunchgrub[testdelays],LegendLabel->"#[delays]"]*)


(* ::Text:: *)
(*Locate position of "r" within degrees("n") being tested*)


(* ::Input:: *)
(*findnprojsintestdegs = (Flatten[Position[crunchgrub[testdegs],liftgrub[nprojs]]])[[1]];*)


(* ::Item:: *)
(*Delays reduce the jump gap from 10^15\[LongRightArrow]10^5 : LTI 1, LTI 2, LTI 3*)


(* ::Input:: *)
(*tplots=MapThread[BoxWhiskerChart[splog10[#1],ChartStyle->#2,ChartLabels-> crunchgrub[testdegs],ImageSize->Large,PlotRange->All,LabelStyle->Directive[Black, Medium]]&,{ensembledat,basicolourlist}];*)
(*(*Manipulate[tplots[[i]],{i,1,25,1}]*)*)


(* ::Subsection:: *)
(* \!\(\*SubsuperscriptBox[\(\[Delta]\), \(rel\), \(\[Mu]\)]\) Vs n*)


(* ::Input:: *)
(*transitplot=(Show[##,ParametricPlot[{findnprojsintestdegs,t},{t,-17,1},PlotStyle->Hue[0.8],PlotRange->All,ImageSize->Large,LabelStyle->Directive[Black,Medium]],PlotRange-> All]&)@@(tplots);*)
(*Show[transitplot,FrameLabel->{{HoldForm[Subscript[Log, 10][\!\(\*SubsuperscriptBox[\(\[Delta]\), \(rel\), \(\[Mu]\)]\)]],None},{HoldForm[n],None}},PlotLabel->None,LabelStyle->{Directive[Black, Medium]}]*)


(* ::Subsubsection:: *)
(*Subscript[\[Sigma], Tail] Vs n*)


(* ::Input:: *)
(*(Show[##,ParametricPlot[{findnprojsintestdegs,t},{t,-17,1},PlotStyle->Hue[0.8],PlotRange->All,ImageSize->Large,LabelStyle->Directive[Black,Medium]],PlotRange-> All]&)@@(MapThread[BoxWhiskerChart[splog10[#1],ChartStyle->#2,ChartLabels-> crunchgrub[testdegs],ImageSize->Large,PlotRange->All,LabelStyle->Directive[Black, Medium]]&,{ensemblechoppeddat,basicolourlist}])*)


(* ::Section:: *)
(*Eigenvalues from multiple time traces : Algorithm 6.1*)


(* ::Subsection:: *)
(*Look at \!\(\*SubsuperscriptBox[\(\[Delta]\), \(rel\), \(\[Mu]\)]\) Vs n and pick the smallest x-coordinate with a large median / the dimension of the underlying system if known*)


(* ::Item:: *)
(*Reduce your guess if indicators below are unsatisfactory*)


(* ::Input:: *)
(*AssociateTo[crunchgrub,rmin-> Input["Please enter your guess of rmin ('r')"](* 7*)];*)


(* ::Subsection:: *)
(*Determining chosendeg*)


(* ::Item:: *)
(*For "n" ranging from cdegmin\[LongRightArrow] cdegmax,  *)


(* ::Subitem:: *)
(*Find the "r" (rmin) estimated eigenvalues using Algorithm 6.1*)


(* ::Subsubitem:: *)
(*True eigenvalues known : Find the Hausdorff distance of the estimated lot to the truth*)


(* ::Subsubitem:: *)
(*If not, compute the pair-wise Hausdorff distances over the range of "n" under purview*)


(* ::Subsubsection:: *)
(* Subscript[Log, 10][Subscript[\[Rho], Hausdorff](\!\(\*OverscriptBox[*)
(*SubscriptBox[\(\[Nu]\), \(n\)], \(^\)]\),\[CapitalLambda])] Vs n : Distance to truth when true eigenvalues are known*)


(* ::Input:: *)
(*getetvdepwrtcdeg[vals,crunchgrub,crunchgrub[rmin]+1 (* 0 would squish the variations with increasing "n"*),crunchgrub[maxn],trajgrub[discevals]]*)


(* ::Subsubsection:: *)
(*Subscript[Subscript[Log, 10][Subscript[\[Rho], Hausdorff](\!\(\*OverscriptBox[*)
(*SubscriptBox[\(\[Nu]\), \(r + i\)], \(^\)]\),\!\(\*OverscriptBox[*)
(*SubscriptBox[\(\[Nu]\), \(r + j\)], \(^\)]\))], i,j>= 1] : Pair-wise Hausdorff distances when true eigenvalues are unknown*)


(* ::Input:: *)
(*getetvdepwrtcdeg[vals,crunchgrub,crunchgrub[rmin]+1 (* 0 would squish the variations with increasing "n"*),crunchgrub[maxn]]*)


(* ::Subsection:: *)
(*Choose chosendeg  and rmin from the earlier diagnostics*)


(* ::Item:: *)
(*chosendeg: There is no set method to choosing this. Following are some suggestions*)


(* ::Subitem:: *)
(*Smallest x-coordinate with low spread in \!\(\*SubsuperscriptBox[\(\[Delta]\), \(rel\), \(\[Mu]\)]\) *)


(* ::Subitem:: *)
(*Easy when you know the true eigenvalues : Pick a point in the first plot that has a low y-coordinate*)


(* ::Subitem:: *)
(*When eigenvalues are unknown (as in the cavity flow), look for a value of 'n' in the second plot that*)


(* ::Subsubitem:: *)
(*Has a small value  in the (1,n)th cell*)


(* ::Subsubitem:: *)
(*Is small for (n,n+1), (n,n+2), (n,n+3).... for as long as possible (i.e.  choosing more estimates of the eigenvalues doesn't change the n-th estimate greatly)*)


(* ::Item:: *)
(*rmin : Update if needed*)


(* ::Input:: *)
(*AssociateTo[crunchgrub,{chosendeg->Input["Please enter your guess of chosendeg ('n' in Algo 6.1)"](*17*),rmin-> Input["Please update your guess of rmin ('r')"](*7*)}];*)
(*{prunedrootqfs,estruevals}=sneakyestimatetruevals[vals,crunchgrub,"verbose"];*)


(* ::Subsubsection:: *)
(*Hausdorff distance between estimated and true eigenvalues*)


(* ::Input:: *)
(*hausdorffdist[estruevals (* Estimated *),trajgrub[discevals] (* True eigenvalues *)]*)


(* ::Subsubsection::Closed:: *)
(*Minutia*)


(* ::Input:: *)
(*oldvals = vals;*)
(*ncases=4;*)


(* ::Subsubsection::Closed:: *)
(*Generate coordinates for heatmaps*)


(* ::Input:: *)
(*plotcoords = Flatten[Outer[{##}&,plotgrub[testdelays],plotgrub[testdegs]],1];*)
(*crunchcoords = Flatten[Outer[{##}&,crunchgrub[testdelays],crunchgrub[testdegs]],1];*)


(* ::Input:: *)
(*vals2plotopsVScases[vals_,crunchcoords_,ncases_]:=Module[{casewiserates,opsVScases,funnyopsVScases,plotopsVScases},*)
(*casewiserates = (* Case, IC,del,deg*)Table[ Map[(#[ratequads])[[All,All,i]]&,vals],{i,ncases}];opsVScases = Table[casewiserates[[j,All,All,All,i]],{i,3}(* PDF, rperror, respercent *),{j,ncases}(* Vanilla, MS, MS + Pres, MS + Delay *)];*)
(*funnyopsVScases = Map[makecratesfunny[Mean,#&,#]&,Rest[opsVScases],{2}];*)
(*plotopsVScases = Map[ldplot[crunchcoords,splog10@#]&,funnyopsVScases ,{2}];*)
(*MatrixForm[plotopsVScases\[Transpose]]*)
(*]*)


(* ::Section:: *)
(*Theorem checks*)


(* ::Item:: *)
(*Each case generates a 4 x 4 table of heat-maps, each plot describing the variation of an averaged quantity with respect to #[delays] and n*)


(* ::Subitem:: *)
(*Columns : Subscript[\[Rho], Subset], Subscript[\[Delta], Trivial]*)


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
(**)
(*vals =ParallelMap[ms1shot4trajvariations[trajgrub,liftgrub,priorsgrub,testpriorsgrub,crunchgrub,#]&,{listoICs,listotseries,oldvals}\[Transpose]];*)
(*vals2plotopsVScases[vals,crunchcoords,ncases]*)


(* ::Subsection:: *)
(*True eigenvalues: \[CapitalLambda]*)


(* ::Input:: *)
(*liftgrub[truevals] = N@trajgrub[discevals];*)
(**)
(*valshonest =ParallelMap[ms1shot4trajvariations[trajgrub,liftgrub,priorsgrub,testpriorsgrub,crunchgrub,#]&,{listoICs,listotseries,oldvals}\[Transpose]];*)
(*vals2plotopsVScases[valshonest,crunchcoords,ncases]*)
