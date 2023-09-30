(* ::Package:: *)

  BeginPackage[ "Atomizer`"]

  gridTree::usage = 
	"gridTree[e] displays a compact tree representation of an expression."
  lispify::usage = "lispify[e] converts a Mathematica expression to a list in Polish Prefix Notation.a"

  Begin["`Private`"]

SetAttributes[dpyNullary,HoldAll];
dpyNullary[ex_] :=
  Grid[{{ex, ""}},
   Frame      -> {All,False},
   Alignment  -> Left,
   Background -> {{LightOrange,{LightYellow}}}];

SetAttributes[dpyMultiary,HoldAll];
dpyMultiary[key_, vals_] :=
   With[{c = Length @ vals},
    Module[{
      spans = Table["", {c}],
      slot  = Floor[(1+c)/2]},
     spans[[slot]] = key;
     Grid[MapThread[List,{spans, vals}],
       Frame      -> {All,False},
       Alignment  -> Left,
       Background -> {{LightOrange,{LightGreen}}}]
      ]];

SetAttributes[dpyAtom,HoldAll];
dpyAtom["String"[e_]]:=
	Grid[{{Style[e,Bold]}},
	Frame     -> All,
	Alignment -> Left,
	Background-> LightYellow];
dpyAtom[ex_] := 
  Grid[{{Style[ex,Bold]}},
   Frame      -> All,
   Alignment  -> Left,
   Background ->
     Switch[Head @ ex,
       String,   LightYellow,
       Symbol,   LightPurple,
       Integer,  LightBlue,
       Real,     LightBlue,
       Rational, LightBlue,
       Complex,  LightBlue,
       _,        Red]];

SetAttributes[lispify,HoldAll];
lispify[c:Complex[_,_]]:=c;
lispify[r:Rational[_,_]]:=r;
lispify[h_[args___]]:=Prepend[lispify/@Unevaluated@{args},lispify[Unevaluated@h]];
lispify[s_(*?AtomQ*)]/;AtomQ[Unevaluated[s]]:=s;(* LEAK!!!  todo *)

SetAttributes[stringulateLisp,HoldAll];
stringulateLisp[l_List]:=stringulateLisp/@l;
stringulateLisp[e_Symbol]:=ToString[Unevaluated[e]];
stringulateLisp[e_String]:="String"[Unevaluated[e]];
stringulateLisp[e_]:=e;(* LEAK!!! todo *)

SetAttributes[gridTree,HoldAll];
gridTree2[l_List]/;Length[l]===1&&AtomQ[First[l]]:=dpyNullary[First[l]];
gridTree2[l_List]/;Length[l]===1:=dpyNullary[gridTree2[First[l]]];
gridTree2[l_List]/;AtomQ[First[l]]:=dpyMultiary[First[l],gridTree2/@Rest[l]];
gridTree2[l_List]:=dpyMultiary[gridTree2[First[l]],gridTree2/@Rest[l]];
gridTree2[e_String]:=dpyAtom[ToExpression[e]];
gridTree2[e_]:=dpyAtom[e];
gridTree2[___]:=Throw["GRIDTREE: CATASTROPHE"];
gridTree[e_]:=gridTree2[stringulateLisp@lispify[e]];

  End[]

  EndPackage[]


