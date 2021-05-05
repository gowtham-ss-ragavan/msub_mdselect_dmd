(* ::Package:: *)

(* ::Title:: *)
(*logistics*)


(* ::Subsubsection::Closed:: *)
(*Update numerical checks using estimatruevals*)


ms1shot4trajvariations[trajgrubIP_, liftgrub_, priorsgrub_, testpriorsgrub_, crunchgrub_, {tinitIP_, rawicsIP_, valselem_}] := Module[{trajgrubOP},
   trajgrubOP = trajgrubIP;
   AssociateTo[trajgrubOP, {tinit -> tinitIP, rawics -> rawicsIP}];
   meansuboneshot[trajgrubOP, liftgrub, priorsgrub, testpriorsgrub, crunchgrub, valselem]
   ];


(* ::Subsubsection::Closed:: *)
(*Generating an IC with specific non-zeros*)


(* ::Item:: *)
(*For data from a finite dimensional Koopman invariant subspace, this can help select the modes to be in the trajectory*)


(* ::Item:: *)
(*Typically, set it to the max possible value to get all modes (as is assumed in theory)*)


getic[ssdim_,nprojs_]:=(*Permute the randomly generated so we spread out the zeros from the end*)RandomSample@(PadRight[(RandomReal[{1,5}(* Trying to keep the ICs within the same order *),nprojs]*RandomChoice[{-1,1}(* Flip signs at random*),nprojs]),ssdim (*Make of correct size*)]);


(* ::Subsubsection::Closed:: *)
(*A plotters digression*)


optchoices2seq[optheads_, optchoices_] := Sequence @@ (Thread[optheads -> optchoices]);


convexcombo[t_,{a_,b_}]:=a + (b-a)t;


hline[y0_,xlimits_]:={convexcombo[#,xlimits],y0}&;
vline[x0_,ylimits_]:={x0,convexcombo[#,ylimits]}&;


makelines[coordlist_,limitlist_,mkrfun_]:=MapThread[mkrfun,{coordlist,limitlist}];


hlines[y0list_,xlimitlist_]:=makelines[y0list,xlimitlist,hline];
vlines[x0list_,ylimitlist_]:=makelines[x0list,ylimitlist,vline];


head4apcore[flavour_ ]:= Switch[flavour,(*Generate function that can MapThread over multiple datasets for the chosen flavour *)
"listplot",ListPlot[Transpose[Most[#1],{3,1,2}],PlotMarkers->#1[[-1]],##2]&,
"parametricplot",ParametricPlot[Evaluate[Through[(#1[[1]])[t]]],{t,0,1},PlotStyle->#1[[2]],##2]&
];


autoplotcore[head_,dataset_,optchoices_]:=Module[{optheads,optsequence},
optheads = {AxesLabel,PlotRange,LabelStyle,ImageSize};
(* Convert choices into a sequence that can be feed into head *)
optsequence = optchoices2seq[optheads,Rest@optchoices];
head[dataset,optsequence]
];


autofigure[style_,grub__]:=autoplotcore[head4apcore[style],grub];


(* ::Item:: *)
(*listplotdat : Typically, data from simulation *)


(* ::Item:: *)
(*parplotdat: Expected values / Lines for interesting values*)


(* ::Item:: *)
(*optchoices: Common settings for  the plot: See autoplotcore for these*)


gildedfigure[listplotdat_,parplotdat_,optchoices_]:=Labeled[Show[MapThread[autofigure[##,optchoices]&,{{"listplot","parametricplot"},(* See head4mpthrd for how to format this data *){listplotdat,parplotdat}}],PlotRange-> All],optchoices[[1]],Top, LabelStyle->optchoices[[3]]];


(* ::Subsubsection::Closed:: *)
(*How to create dummy vals for the sake of testing*)


getdummycvec[truevals_,deg_]:=Module[{locvals,n,randrange = {-1-I,1+I},seedvals},
n = Length[truevals];
seedvals = Join[truevals,RandomComplex[randrange,Abs[n-deg]]];
locvals = If[deg<n,(* Chose a random subset if you are smaller than n *)RandomSample[seedvals,deg],(* Contains truevals and is of length deg*)seedvals];
(* Snip the -1 at the end *)Most[(* Negative to get in Companion matrix format*)-1*roots2coeffs[locvals]]
];


getdummyvalselement[delays_,degs_,truevals_]:=<|cquads-> Outer[(* Braces essential as this is supposed to be a set of cvectors from all 4 cases *){getdummycvec[truevals,#2]}&,delays,degs],equads-> Outer[(* Braces essential as this is supposed to be a set of cvectors from all 4 cases *){{1,1,1}}&,delays,degs]|>;


getdummyvals[delays_,degs_,truevals_,nics_]:=Table[(* Generate Association for each IC *)Unevaluated[getdummyvalselement[delays,degs,truevals]],{i,nics}];


(* ::Subsection::Closed:: *)
(*Extracting parts of vals (Data structure used in scripts)*)


(* ::Subsubsection::Closed:: *)
(*c* for plain DMD @ a particular degree*)


vals2vanillactes[vals_,degindex_]:=Flatten[Map[(#[equads])[[All,degindex,1]]&,vals],1];


(* ::Subsubsection::Closed:: *)
(*Residue info for plain DMD @ a particular degree*)


vals2vanillacvecs[vals_,degindex_]:=Flatten[Map[(#[cquads])[[All,degindex,1]]&,vals],1];


(* ::Subsubsection::Closed:: *)
(*Prunes the delays so that the rest have a consistent jump index*)


(* ::Item:: *)
(*Needs crunchgrub to have the info on the minimum number of delays to take in delaymin*)


keepthemgoodelays[vals_, crunchgrub_] := Module[{goodstart},
   (* Find the index within testdelays that corresponds to delaymin *)
   goodstart = (Position[crunchgrub[testdelays],crunchgrub[delaymin]])[[1, 1]];
   vals[[All(*ICs*), All(*Keys within each Association *),(*Delays*)goodstart ;; -1]]
   ];
