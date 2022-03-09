(* ::Package:: *)

(*
 Wolfram Language scrapbook file.

 The contents can be evaluated by selecting the 
 file, right-clicking, and choosing Wolfram/Run
 or Wolfram/Debug.
*)

genericpickysvd[x_, sigpicker_] :=
    Module[{usv, rawsigmas, sigmas, condandterror},
        (*Get SVD of x*)
        (*UpTo: To compute all Subscript[\[Sigma],i] so that we can threshold ourselves*)
        usv = SingularValueDecomposition[x, UpTo@Min[Dimensions[x]]];
        (*Here we choose the singular values that are considered and those that are neglected*)
        rawsigmas = Diagonal[usv[[2]]];
        sigmas = sigpicker@rawsigmas;
        (*Find unthresholded and thresholded \[Kappa](x) and the largest Subscript[\[Sigma],i] not included*)
        condandterror = getcondandtail[rawsigmas, Length@sigmas];
        (*Choose the corresponding rows of U^H and V^H*)
        {thinsvd[usv, sigmas], condandterror}
    ];

pickysvd[x_, sigtols_] := 
  genericpickysvd[x, sigthresholder[#, sigtols] &];

lowrankapprox[x_, sigpicker_] := 
  usvdot@((genericpickysvd[x, sigpicker])[[1]]);

(* Removed tinit -  Carries noise too *)

ms1shot4trajvariations[trajgrubIP_, liftgrub_,  crunchgrub_, {rawicsIP_,noiseIP_, valselem_}] := Module[{trajgrubOP},
   trajgrubOP = trajgrubIP;
   AssociateTo[trajgrubOP, { rawics -> rawicsIP, opnoise -> noiseIP }];
   meansuboneshot[trajgrubOP, liftgrub, crunchgrub, valselem]
   ];
  
 

getcquad[data_,\[Lambda]_,sigtols_,restols_,truevals_,nfunda_,crows_,ndelays_,deg_,meff_, noisyQ_]:=Module[{datapartsflavs},
(* First create the various inputs for each flavour of Companion DMD *)
datapartsflavs = prep4cases[data,\[Lambda],truevals,sigtols, restols,nfunda, ndelays, deg, meff, noisyQ];
MapThread[getcstar[#1,\[Lambda],sigtols,restols,#2,#3,crows,ndelays]&,datapartsflavs]\[Transpose](* [[Attributes, Cases ]]*)
];

prep4cases[data_, \[Lambda]_, truevals_, sigtols_, restols_, nfunda_, ndelays_, deg_, meff_, noisyQ_] :=
    Module[ {zmat, zlambdams, zlambdamsdel, n, evalparts, flavours, dataorgs,m,msdelgrub,msdelgrubsansmean},
    	(*Model parameters *)
        n = deg;
        (* True length of observables*)
        m = Length[data[[1]]]/(ndelays + 1);
        (*Actual data matrices to be parsed*)
        zmat = Take[Transpose[data[[1]]], -(n + 1)] // Transpose;
        zlambdams = getmusubeddatamat[zmat, \[Lambda]];
        (* Need to remove the mean, and then take time delays *)
        (* Pick out the time series as it was before the delays, of appropriate length *)
        msdelgrub = Take[Transpose[ data[[1,Range[m],All]] ], -(n+1+ndelays)]//Transpose;
        (* Remove temporal mean *)
        msdelgrubsansmean = getmusubeddatamat[msdelgrub, \[Lambda]];
        (* Do delays *)
        zlambdamsdel = tseries2bigdatamat[(msdelgrubsansmean)\[Transpose], ndelays - 1];
        evalparts = getevalpartitions[truevals, nfunda, n, \[Lambda], restols];
        (* Update the eigenvalue target, since the mean was using a matrix of length n+ndelays+1 *)
        (* True without rate always *)
        evalparts[[-2]] = (getevalpartitions[truevals, 0, n+ndelays, \[Lambda], restols])[[-1]];
        (* True with rate always *)
        evalparts[[-1]] = (getevalpartitions[truevals, n+ndelays+1, n+ndelays, \[Lambda], restols])[[-1]];
        flavours = {"vanilla", "ms", "msdelay", "msdelay"};
        dataorgs = {zmat (*Usual*), zlambdams (*Mean subtracted*), 
          zlambdamsdel (*Remove mean-Preserve \[Lambda]*), 
          zlambdamsdel(*Remove mean,then delay*)};
        {dataorgs, evalparts, flavours}
    ] /; (! noisyQ)




prep4cases[{data_, noise_}, \[Lambda]_, truevals_, sigtols_, restols_,
    nfunda_, ndelays_, deg_, meff_, noisyQ_] :=
    Module[ {filterNgrub, zpure, znoisy, zTLS, zfilter, znoiseresist, n, evalparts, flavours, dataorgs, filter = <||>, grub = <||>, sigpicker},
        n = deg;
        (*Take the first meff and last n+1 columns please*)
        filterNgrub = {Take[# // Transpose, meff] // Transpose, 
           Take[# // Transpose, -n - 1] // Transpose} &;
        {filter["data"], grub["data"]} = filterNgrub[data];
        {filter["noise"], grub["noise"]} = filterNgrub[noise];
        (*Split data and noise into filter and actual grub*)
        (*Companion: Pure data*)
        zpure = grub["data"];
        (*Noisy companion:Data+Noise*)
        znoisy = zpure + grub["noise"];
        (*TLS Companion: Data+Noise-->SVD thresholding*)
        sigpicker = 
         If[ Length[truevals] > 0,
         	(*Retain only as many singular values as there are Koopman EFs*)
         	(*Fold[sigthresholder, {#, UpTo[Length[truevals]], sigtols}] &*)
             sigthresholder[#, UpTo[Length[truevals]]] &,
             sigthresholder[#, sigtols] &
         ];
        zTLS = lowrankapprox[znoisy, sigpicker];
        (*Noise resistant DMD: Noise added to both data matrices --> Inner product with the average in front*)
        zfilter = filter["data"] + filter["noise"];
        znoiseresist = (1/(ndelays + 1))*
          ConjugateTranspose[zfilter] . znoisy;
        evalparts = ConstantArray[{truevals, {}}, 4];
        flavours = ConstantArray["ms", 4];
        dataorgs = {zpure, znoisy, zTLS, znoiseresist};
        {dataorgs, evalparts, flavours}
    ];


updatednumchks[data_,\[Lambda]_,truevals_,restols_,sigtols_,nfunda_,ndelays_, deg_,meff_, noisyQ_,locquads_]:=
Module[{dataorgs,evalparts,flavours},
{dataorgs,evalparts,flavours} = prep4cases[data,\[Lambda],truevals,sigtols,restols,nfunda, ndelays, deg, meff, noisyQ];
MapThread[numchks4msubthrms[zmat2xmat[#1],#2,#3,sigtols]&,{dataorgs,evalparts,locquads}]
];




meansuboneshot[grubsNprior__Association]:=Module[{grubs,prior,args},
{prior,grubs} = TakeDrop[{grubsNprior},-1];
args = Join[{maxdelays,rawics,opnoise,noisyQ,meff,sigtols,restols,cmatseed,crows,crank,testdelays,testdegs,rate2sub,nfunda,truevals}/.(Normal[Join@@grubs]),prior];
meansuboneshot@@args
];



(*------Generic meansuboneshot that should work for both deterministic and stochastic cases ---------*)

meansuboneshot[ndelays_ ,rawics_,opnoise_,noisyQ_, meff_, sigtols_,restols_,cmatseed_,crows_,crank_,testdelays_,testdegs_,rate2sub_,nfunda_,truevals_,ipgrub_]:= 
Module[{obsdim,datamat,noisemat,tseries,cmat,opgrub,crequads,crefun,matcarvefun,deldegmat,newchks,newrqs},
(*-----------------------------*)
(*------PRELIMINARIES-------*)
(*-----------------------------*)

(* Parameters that we'll need later on*)
(* obsdim = r in the manuscript *)
obsdim = Length[rawics[[1]]]; 


(* cmat ~ m x r = crows x obsdim *)

cmat = If[cmatseed == {}(* If unspecified, create using getrandmat *), 
	usvdot@((* randmat's OP is in usv form*)getrandmat[crows,obsdim,{crank,{}}]),
	cmatseed (* Use given Cmatrix*)];



(* Create the ICs*)
tseries = rawics . Transpose[cmat] ;



(* The maximally delayed Z matrix *)

{datamat, noisemat} = Map[tseries2bigdatamat[#,ndelays-1]&, {tseries, opnoise}];  





(* Compute c-r-e quads over the desired set of delays, trajectory length / system order *)


matcarvefun = basicmatcarve[{datamat,noisemat}, crows, #1, #2,meff]&;

If[ipgrub==<||>,
(*------------- First computational pass : No priors ----------------*)




crefun =getcquad[matcarvefun[#1,#2],rate2sub,sigtols,restols,truevals,nfunda,crows,#1,#2, meff,noisyQ]&;








crequads=Transpose[Outer[crefun,testdelays,testdegs],{2,3,1,4}];
opgrub=<|cquads-> crequads[[1]],ratequads-> crequads[[2]],equads-> crequads[[3]] |>,


(*------------- Updating numerical checks with new eigenvalue estimates ----------------*)


opgrub = ipgrub;
(* Compute rperror and respercent with truevals whiich is likely a product of estimatruevals*) 

crefun = updatednumchks[matcarvefun@@#1,rate2sub,truevals,restols,sigtols,nfunda,#1[[1]],#1[[2]], meff, noisyQ, #2]&;


deldegmat = Outer[{##}&,testdelays,testdegs];
newchks = MapThread[crefun,{deldegmat,ipgrub[cquads]},2];
(*Update opgrub with newchks appropriately *)
newrqs = opgrub[ratequads];
newrqs[[(* Delays *)All,(* Degrees *)All,(* Cases *)All,{2,3}]]=newchks;
opgrub[ratequads]=newrqs
];





opgrub
];


keepthemgoodelays[vals_, crunchgrub_] := Module[{goodstarts},
   (* Find the index within testdelays that corresponds to delaymin *)
   goodstarts = Flatten@Position[Thread[crunchgrub[testdelays] >= crunchgrub[delaymin] ],True];
   vals[[All(*ICs*), All(*Keys within each Association *),(*Delays*)goodstarts]]
   ];


(*Intend to create a function that will save variables with plots as their members

SICK of those annoying paint operation to merge these things

Could potentially use this in the future.

I/P: Variable, Function that names the files when given the index, file format

*)

savetheseplots[plotvar_,plotnamer_,format_]:= MapIndexed[Export[StringJoin[plotnamer[#2[[1]]],".",format],#1]&,plotvar];

nametheseplots[i_,casestrings_]:=Module[{systemname = crunchgrub[savefile],nameroot,casename},
nameroot = FileNameJoin[{crunchgrub[plotdir],systemname}];
casename = casestrings[[i]];
StringJoin[nameroot,"_",casename]
];


sparsifyChartLabels[labels_, spacing_] :=
  MapIndexed[
   Switch[
     Mod[#2[[1]], spacing],
     1, #1,
     _, " "
     ] &, labels];
     




(*1 : Good, 0: Bad*)

kmdQuality[{valswts_, rperror_, respercent_}] := 
  1 - Max[1 - 10^-rperror, respercent];
evalsQuality[{valswts_, rperror_, respercent_}] := 10^-rperror;
explainQuality[{valswts_, rperror_, respercent_}] := 1 - respercent;

 
prefixstrlist[prefix_,strlist_]:=Map[StringJoin[prefix,#]&,strlist];

(* Process you data before plotting - Makes labelling more intuitive *)


basicBWplot[allXlabels_, ensembledat_, 
   colourlist_, {vertmin_, vertmax_}, vertcolour_] := 
  Module[{tplots, xlabelspacing = 3, sparseXlabels},
   (* Plot only a subset of xlabels to avoid cluttering for larger \
fonts *)
   (*Replace those whose index %xlabelspacing != 1 with None *)
   sparseXlabels = sparsifyChartLabels[allXlabels, xlabelspacing];
   tplots = MapThread[
     BoxWhiskerChart[#1, {{"Whiskers", Transparent}, {"Fences", 
         Transparent}, {"MedianMarker", "o", #2}}, ChartStyle -> #2, 
       ChartLabels -> sparseXlabels, ImageSize -> Large, 
       PlotRange -> All, LabelStyle -> Directive[Black, Medium], 
       Joined -> True, BarSpacing -> 0.6 (* 
       Allows the Join of Medians to be discerned, 
       if reader is interested *),  PlotTheme -> "Scientific"] &,
     {ensembledat, colourlist}];
   ( 
     Show[##, 
       ParametricPlot[{findnprojsintestdegs, t}, {t, vertmin, 
         vertmax}, PlotStyle -> Directive[vertcolour, Dashed](*Hue[
        0.8]*), PlotRange -> All, ImageSize -> Large, 
        LabelStyle -> Directive[Black, Medium]], PlotRange -> All] &
     ) @@ (tplots)
   ];

stdBWplot[allXlabels_, ensembledat_, colourlist_, colorlegend_, 
      xlabel_, ylabel_, vertbounds_, vertcolour_, placement_] :=
    Module[ {transitplot},
        transitplot = 
          basicBWplot[allXlabels, ensembledat, colourlist, vertbounds, 
            vertcolour];
        (*Put the correct legend*)
        (Legended[#, 
        Placed[colorlegend, (* Nice delay plots are better displayed if the legend were within the SE corner *)
        placement]] &)@(*On this bunch of plots that have been \
        overlaid*)Show[transitplot, 
        FrameLabel -> {{ylabel, None}, {xlabel, None}}, 
        PlotLabel -> None, LabelStyle -> {Directive[Black, 20]}]
    ];
   
stdBWplot[allXlabels_, ensembledat_, colourlist_, colorlegend_, 
      xlabel_, ylabel_, vertbounds_, vertcolour_] :=
    stdBWplot[allXlabels, ensembledat, colourlist, colorlegend, 
        xlabel, ylabel, vertbounds, vertcolour, After];

kmdplots[allXlabels_, key_, fun_, vals_, funname_, colourlist_(**), 
   colorlegend_(**), vertbounds_, vertcolour_] := Module[{funnyvals, plots},
      (* Cases, Del, Deg, ICs *)
      funnyvals = Transpose[
          Map[
            Map[fun, #[key], {3}] &,
            vals],
          4 <-> 1
          ];
      plots = 
        Map[stdBWplot[allXlabels, #, colourlist , colorlegend, "Companion-order", 
              funname, vertbounds, vertcolour] &, funnyvals];
      plots
      ]; 

kmdplots[allXlabels_, key_, fun_, vals_, funname_, colourlist_(**), 
   colorlegend_(**)] := 
   kmdplots[allXlabels, key, fun, vals, funname, colourlist, 
   colorlegend(**), {0,1}, Hue[0.8]];
   


