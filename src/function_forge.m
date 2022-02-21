(*
 Wolfram Language scrapbook file.

 The contents can be evaluated by selecting the 
 file, right-clicking, and choosing Wolfram/Run
 or Wolfram/Debug.
*)

genericpickysvd[x_, sigpicker_] :=
    Module[ {usv, rawsigmas, sigmas, condandterror},
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
    Module[ {zmat, zlambdams, zlambdamsdel, n, evalparts, flavours, dataorgs},
        n = deg;
        zmat = Take[Transpose[data[[1]]], -(n + 1)] // Transpose;
        zlambdams = getmusubeddatamat[zmat, \[Lambda]];
        zlambdamsdel = 
         Join @@ Map[#[zlambdams\[Transpose]]\[Transpose] &, {Most, Rest}];
        evalparts = 
         getevalpartitions[truevals, nfunda, n, \[Lambda], restols];
        flavours = {"vanilla", "ms", "mspres", "msdelay"};
        dataorgs = {zmat (*Usual*), zlambdams (*Mean subtracted*), 
          zlambdams (*Remove mean-Preserve \[Lambda]*), 
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
             Fold[sigthresholder, {#, UpTo[Length[truevals]], sigtols}] &,
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
tseries = rawics.Transpose[cmat] ;



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

