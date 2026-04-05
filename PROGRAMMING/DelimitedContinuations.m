(* ::Package:: *)

(* ::Input:: *)
(*ClearAll[reset,shift];*)
(*reset[h_[As___, shift[k_, E_], Bs___]]:=*)
(*Block[{K=x\[Function]reset[h[As,x,Bs]]},*)
(*Print[<|"As"->{As},"h"->h,"k"->k,"E"->E,"Bs"->{Bs}|>];*)
(*reset[E/.{k->K}]];*)
(*reset[E_]:=E;*)


(* ::Text:: *)
(*The following shows that reset is transparent to expressions that don't contain shifts.*)


(* ::Input:: *)
(*reset[{1,2,3}]*)


(* ::Text:: *)
(*The following shows that the pattern recurses properly over multiple shifts in the argument of reset.*)


(* ::Input:: *)
(*reset[{1,shift[k1,k1[E1]],2,shift[k2,k2[E2]],3}]*)


(* ::Text:: *)
(*The following shows that reset's contents are arbitrary, and can still be rewritten:*)


(* ::Input:: *)
(*reset[arbitrary[1,shift[k1,k1[E1]],2,shift[k2,k2[E2]],3]]*)


(* ::Text:: *)
(*This shows calling the continuation in the body of the shift.*)


(* ::Input:: *)
(*2*reset[1+shift[k,k[5]]]*)


(* ::Text:: *)
(*This shows calling the continuation multiple times in the body of the shift.*)


(* ::Input:: *)
(*reset[2*shift[k,k[k[4]]]]*)


(* ::Text:: *)
(*Here are a couple of examples from https://goo.gl/Mf8K0y. *)


(* ::Input:: *)
(*reset[{shift[k, Prepend[k[{}],1]]}]*)


(* ::Input:: *)
(*reset[{*)
(*shift[k, Prepend[k[{}],1]],shift[k, Prepend[k[{}],2]]}]*)


(* ::Input:: *)
(*ClearAll[yield];*)
(*yield[x_]:=shift[k,Prepend[k[{}],x]];*)


(* ::Input:: *)
(*reset[{yield[1],yield[2],yield[3]}]*)



