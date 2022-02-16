(* ::Package:: *)

(* ::Title:: *)
(*logistics*)


(* ::Subsubsection::Closed:: *)
(*Process a time series before superbasislift*)


tseries2bigdatamat[tseries_,ndelays_]:=Map[Flatten,Partition[tseries,ndelays+2,1]]\[Transpose];
(* You should have ndelays + 1: An additional 1 is required to compensate for the loss in getXtendedmats, via liftedY *)


(* ::Subsubsection::Closed:: *)
(*{x_n + I y_n} -->  {{x_n, y_n}}*)


complexto2D[cvaluelist_] := Through[{Re, Im}[cvaluelist]]\[Transpose];


(* ::Subsubsection:: *)
(*Complex numbers colored by some attribute*)


(* ::Item:: *)
(*att2pow:  That which turns attributes into powers ~ Positive numbers whose magnitude correlates with importance*)


powplaneplot[evalattpairs_, att2pow_] := Module[{evals, atts, pows, powbounds, huenums, lpgrub},
   {evals, atts} = evalattpairs\[Transpose];
   pows = Map[att2pow, atts];
   powbounds = Through[{Min, Max}[pows]];
   huenums = Rescale[pows, powbounds, {0, 0.7}(* Cycle from Red to Blue as the power increases *)];
   lpgrub = MapThread[Style[#1, Hue[#2]] &, {complexto2D@evals, huenums}];
   ListPlot[lpgrub, PlotStyle -> PointSize[Medium], ImageSize -> Large]
   ];


(* ::Subsubsection::Closed:: *)
(*Z -> X*)


(* ::Text:: *)
(*Select all but the last snapshot in a data matrix*)


zmat2xmat[zmat_]:=zmat[[All,1;;-2]];


(* ::Subsubsection::Closed:: *)
(*Pick that returns Nothing instead of {}*)


politePick[list_,ttable_]:=Module[{buffer},
buffer = Pick[list,ttable];
If[buffer=={},Nothing,buffer]
];


(* ::Subsubsection::Closed:: *)
(*Append a row of -1*)


dresscmat[cmat_]:=ArrayPad[cmat,{{0,1},{0,0}},-1];


(* ::Subsubsection:: *)
(*Rows of a Vandermonde*)


(* Generate a vector of monomials upto power n at \[Lambda] *)
getgpvec[\[Lambda]_,n_]:=Table[splpower[\[Lambda],i],{i,0,n}];


(* ::Subsection::Closed:: *)
(*Comparing two numeric quantities*)


(* ::Subsubsection::Closed:: *)
(*Is this close to 0 ?*)


smallQ[x_]:=TrueQ[Abs[x]<10^-8];


(* ::Subsubsection::Closed:: *)
(*Check two scalars/vectors/matrices and O/P the relative and absolute errors*)


nums4absrelcheck[x1_,x2_]:=Module[{aerror,rerror},
aerror = Norm[Flatten[x2-x1]];
rerror = aerror/Norm[Flatten[x1]];
{rerror,aerror}
];


greaterabsrelcheck[x1_,x2_,{reltol_,abstol_}]:=Module[{aerror,rerror},
(* Note the order of O/P is the same as the order of tolerances - R.A.*)
{rerror,aerror} = nums4absrelcheck[x1,x2];
{TrueQ[aerror > abstol && rerror > reltol],{rerror,aerror}}
];


(* ::Text:: *)
(*Ware this version that uses a 4-component restol*)


greaterabsrelcheck[x1_,x2_,{resrel_,res4lamb_,res4trim_,res4grade_}]:=greaterabsrelcheck[x1,x2,{resrel,res4grade}];


(* ::Subsubsection::Closed:: *)
(*Hausdorff distance between two discrete sets*)


(* ::Item:: *)
(*Cannot handle edge case when one is an empty set *)


hausdorffdist[set1_,set2_](* Compute Hausdorff distance between  2 given sets *):=Module[{dtmat},
dtmat = Outer[Norm[#1-#2]&,N@set1,N@set2,1];
Max@Map[(Max@Map[Min,#])&,{dtmat,dtmat\[Transpose]}]
];


(* ::Subsubsection::Closed:: *)
(*Numeric check for set containment*)


(* ::Item:: *)
(*\[Rho]_{Subset} in draft*)


(* ::Subitem:: *)
(*Modification of Hausdorff distance*)


(* ::Subitem:: *)
(*Algorithm requires both sets be finite*)


pert4lhsinrhs[set1_,set2_](* What is the smallest bubble that must be drawn around each element of set2 to contain set1 : Non-negative and 0 \[DoubleLongLeftRightArrow] set1 \[Subset] set2 *):=Module[{dtmat},
(* Compute the discrete metric matrix *)
dtmat = Outer[Norm[#1-#2]&,N@set1,N@set2,1];
Max@Map[Min,dtmat]
];


(* ::Subsection::Closed:: *)
(*Generating random matrices*)


(* ::Subsubsection::Closed:: *)
(*Generate an UpperTriangular matrix of maximal rank by setting the diagonal entries to 1*)


getuppertriangularmat[m_,n_]:=Module[{rmat,diag,d,correction,randrange = {-10,10}},
rmat = RandomReal[randrange,{m,n}];
diag = Diagonal[rmat];
d = Length@diag;
correction = ArrayPad[DiagonalMatrix[-diag+1],{{0,m-d},{0,n-d}}];
UpperTriangularize[rmat + correction]
];


(* ::Subsubsection::Closed:: *)
(*Get a lower triangular matrix of maximal rank by setting the diagonal to 1*)


getlowertriangularmat[m_,n_]:=getuppertriangularmat[n,m]\[Transpose];


(* ::Subsubsection::Closed:: *)
(*Generate a full rank matrix by explicitly setting up it's LU decomposition*)


getfullrankmat[m_,n_]:=getlowertriangularmat[m,m].getuppertriangularmat[m,n];


(* ::Subsubsection::Closed:: *)
(*Generate a random matrix, with specific ( or #) singular values*)


(* ::Item:: *)
(*O/P in SVD format*)


getrandmat[m_,n_,{r_,sigmas_}]:=Module[{buffer,sigs},
buffer = SingularValueDecomposition[getfullrankmat[m,n]];
sigs = Switch[Length@sigmas,
0,Take[Diagonal[buffer[[2]]],r],
_,sigmas];
thinsvd[buffer,sigs]
];


(* ::Subsubsection::Closed:: *)
(*Converts the O/P of SingularValueDecomposition into the I/P matrix (Multiply them all !!)*)


usvdot[{u_,s_,v_}]:=Dot[u,s,v\[ConjugateTranspose]];


(* ::Subsubsection::Closed:: *)
(*Generating a matrix with specified eigenvalues*)


diagreplace[mat_,correction_]:=mat + DiagonalMatrix[correction - Diagonal[mat]];


eval2amat[evals_]:=Module[{n,rmat,qmat},
n = Length@evals;
rmat = getuppertriangularmat[n,n];
qmat = (SchurDecomposition[rmat\[ConjugateTranspose].rmat])[[1]];
rmat = diagreplace[rmat,evals];
Dot[qmat,rmat,qmat\[ConjugateTranspose]]
];


(* ::Subsection::Closed:: *)
(*Sampling*)


(* ::Subsubsection::Closed:: *)
(*Sampling a trajectory to digitize a continuous curve*)


(* ::Item:: *)
(*putty : Use to convert the trajectory to a different coordinate system*)


(* ::Subitem:: *)
(*Specifically meant for the limit cycle examples*)


(* ::Subitem:: *)
(*Usual systems will use Identity to do this*)


samplethistraj[thead_,{tinit_,tsamp_,ndelays_},putty_]:=Table[putty[thead[tinit + i*tsamp]],{i,0,ndelays+1}];


(* ::Subsubsection::Closed:: *)
(*Generating random ICs*)


(* ::Item:: *)
(*nseeds : Number of points you want*)


(* ::Item:: *)
(*ranges: Specify the range that the values can take in each coordinate*)


getrandpoints[{nseeds_,ranges_}]:=Map[RandomReal[#,nseeds]&,ranges]\[Transpose];


(* ::Subsubsection::Closed:: *)
(*Log_{10},with a cautious positive offset*)


splog10[data_]:=Log10[ data + $MachineEpsilon];
