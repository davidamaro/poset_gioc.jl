#!/home/david/docs/bin/Executables/wolframscript
(* ::Package:: *)

Quiet[<<Combinatorica`];
fileName = $ScriptCommandLine[[2]]
ll = StringLength[fileName]
exportName = StringTake[fileName, ll - 3]<>"png"
data = Import[fileName];
nn = Length[data];
xx = data // AdjacencyGraph // TransitiveReductionGraph // AdjacencyMatrix // Normal;
g = FromAdjacencyMatrix[xx, Type -> Directed];
h = MemoryConstrained[ HasseDiagram[SetVertexLabels[g, Range[1, nn]]], 17597961216];
fig = ShowGraph[h, EdgeDirection -> True, BaseStyle -> {FontSize->20}, VertexColor -> Red, VertexLabelColor -> Black, EdgeColor -> Gray, ImageSize -> 700];
Export[exportName, fig]
Print[exportName]
