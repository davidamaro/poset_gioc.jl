#!/home/david/docs/bin/Executables/wolframscript
(* ::Package:: *)

(*
Funciones de ayuda: comienzo
*)
primer[ll_] := 
 StringReplace[ll, {"University" -> "", "of" -> "", "," -> "", "The" -> ""}]
soloMayus = StringReplace[t_?LowerCaseQ :> ""]
final[st_] := With[{largo = primer[st]},
  If[StringLength[largo] > 10, soloMayus[largo], largo]
  ]
(*
Funciones de ayuda: final
*)

Quiet[<<Combinatorica`];
fileName = $ScriptCommandLine[[2]]
Print["size: ", Length@$ScriptCommandLine]
If[Length@$ScriptCommandLine > 2,
names = final/@Import[$ScriptCommandLine[[3]], "List"];
Print[names];
]
ll = StringLength[fileName]
exportName = StringTake[fileName, ll - 3]<>"png"
data = Import[fileName];
nn = Length[data];
xx = data // AdjacencyGraph // TransitiveReductionGraph // AdjacencyMatrix // Normal;
g = FromAdjacencyMatrix[xx, Type -> Directed];

If[Length@$ScriptCommandLine > 2,
h = MemoryConstrained[ HasseDiagram[SetVertexLabels[g, names]], 17597961216],
h = MemoryConstrained[ HasseDiagram[SetVertexLabels[g, Range[1, nn]]], 17597961216]
];

fig = ShowGraph[h, EdgeDirection -> True, BaseStyle -> {FontSize->20}, VertexColor -> Red, VertexLabelColor -> Black, EdgeColor -> Gray, ImageSize -> 700];
Export[exportName, fig]
Print[exportName]
