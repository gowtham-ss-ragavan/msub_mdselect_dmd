(* ::Package:: *)

(* ::Subsubsection:: *)
(*Creating mein src file*)


Module[{address4source,srcfilenames},
(* Assemble the address of the source directory *)address4source =FileNameJoin[{(* Go one step back *)ParentDirectory@DirectoryName[NotebookFileName[]],(* Choose the src directory *)"src"}];
(* Pick the files that have a particular phrase in them *)
srcfilenames=Flatten@Map[FileNames[__~~#~~__,address4source]&,{"core","logistics","untested"}];
(*Run each .m file listed*)
Scan[Get[#]&,srcfilenames]
];
(* Wipe the temporary variables : Clean freak *)
Remove[address4source,srcfilenames];
