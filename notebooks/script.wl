#!/home/david/docs/bin/Executables/wolframscript
(* ::Package:: *)

fileName = $ScriptCommandLine[[2]]
ll = StringLength[fileName]
exportName = StringTake[fileName, ll - 3]<>"png"
data = Import[fileName]
options = 
  Sequence[VertexStyle -> LightYellow, VertexSize -> 0.7, 
   VertexLabels -> Placed["Name", Center], 
   VertexLabelStyle -> 
    Directive[20, Red, Bold, 
     Italic],(*GraphLayout\[Rule]"CircularEmbedding",*)
   ImageSize -> 1050, EdgeStyle -> Blue];
fig = AdjacencyGraph[data, options] // TransitiveReductionGraph
Export[exportName, fig]
Print[exportName]
