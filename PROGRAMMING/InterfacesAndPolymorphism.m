(* ::Package:: *)

(* ::Title:: *)
(*Polymorphism and Interfaces in Mathematica*)


(* ::Subtitle:: *)
(*Brian Beckman*)
(*7 December 2013*)


(* ::Section:: *)
(*Introduction*)


(* ::Text:: *)
(*Imagine the classic C# Enumerable API, leaving out Reset, and let's implement it for Mathematica Lists and Mathematica hash tables (which are lists of replacement rules subjected to Dispatch). This serves as an illustration of a way to implement interfaces that are polymorphic, that is, can be implemented differently by different concrete providers.*)


(* ::Item:: *)
(*REFERENCE: http://bit.ly/1jxWqHe*)


(* ::Section:: *)
(*Object-Oriented Programming (OOP)*)


(* ::Text:: *)
(*Simulate oop's objects as explicit lists of rules that transform patterns into expressions. When an object goes out of scope, its list of rules is automatically garbage-collected. This implementation avoids stocking the global symbol table with rules for each instance, avoiding explicit reference counting and manual deletion. *)


(* ::Text:: *)
(*Every instance of a class must have a rule for every one of the class's members. Likewise, an object must have a rule for every member of every interface that the object implements. We say that an object provides implementations for its class's and interface's members, or that the object is a provider of its class's and interface's members. In both cases, the object provides implementations by providing rewrite rules. *)


(* ::Subsection:: *)
(*Type-Checking*)


(* ::Text:: *)
(*We insist that an object must be a list of rewrite rules explicitly constructed in a Dispatch table, i.e., an expression with Head Dispatch. Since Mathematica rewrites dispatches back to lists at its own discretion, we cannot be sure that every object is a Dispatch even if we insist that the user explicitly construct objects by invoking Dispatch. Instead, we can only check that objects are either a List or a Dispatch, and that's what the helper function ObjectQ does. This helper is invoked as part of the pattern-matching at call sites of functions that require objects. This invocation at the call site implements run-time type checking.*)


(* ::Input:: *)
(*ClearAll[ObjectQ];*)
(*ObjectQ[this_]:=*)
(*Switch[Head@this,*)
(*Dispatch,True,*)
(*List,(* Ensure every element is a rule *)this===Select[this,With[{k=Head@#},k===Rule||k===RuleDelayed]&],_,False]*)


(* ::Text:: *)
(*To type-check that an object provides a certain class or interface, just check that its list-of-rules or its Dispatch contains rules for each of the required members. Represent types -- either classes or interfaces -- as lists of the patterns -- left-hand-sides -- required for each rule of some object. For example, every provider of IEnumerator must provide a rule with a pattern MoveNext[] -- a void-to-Boolean method -- and a pattern for Current -- a property producing a value. To check that an object provides such rules, just extract the First -- the pattern on the left-hand side -- from each rule of the object, and check that the required patterns are a subset of the result. A subtlety is that actual parameter names must be stripped from type-specs extracted from objects, leaving only the types of parameters. For example, we must convert MoveNext[i_Integer] to MoveNext[_Integer]. We accomplish this by a simple rewrite of Pattern expressions to Blank expressions (though see this known bug in the patterm-matcher: http://bit.ly/1bosOYv).*)


(* ::Input:: *)
(*ClearAll[ProvidesTypeQ,SubsetQ,GetPatterns,StripName];*)
(*SubsetQ[A_,B_]:=Complement[A,B]==={};*)
(*StripName[*)
(*Verbatim[Pattern][nym_,typeSpec:Verbatim[Blank][type___]]]:=typeSpec;*)
(*StripName[hd_[args___]]:=hd@@StripName/@{args};*)
(*StripName[else_]:=else*)
(*GetPatterns[this_?ObjectQ]:=With[{h=Head@this},*)
(*If[h===Dispatch,*)
(*StripName/@First/@this[[1]],*)
(*If[h===List,*)
(*StripName/@First/@this,*)
(*Throw["InvalidOperationException"]]]];*)
(*ProvidesTypeQ[this_?ObjectQ,type_List]:=SubsetQ[type,GetPatterns[this]]*)


(* ::Subsection:: *)
(*Inheritance and Polymorphism*)


(* ::Text:: *)
(*Particular rules could simulate inheritance and polymorphism either by chasing prototype chains (upward inheritance) or dispatching at run time (downward inheritance). We don't need to take a cosmic position on this; each class and interface can do it in its own way. *)


(* ::Subsection:: *)
(*Using Objects*)


(* ::Text:: *)
(*To produce application-level results, apply objects' rules to other, arbitrarily rich work-expressions. This is "rewrite at run time." Contrast with the more familiar "rewrite at compile time," which is how standard oop systems transform work-expressions into actionable, application-level work-code.*)


(* ::Subsection:: *)
(*Overload Resolution; Member Lookup*)


(* ::Text:: *)
(*In the current design, member lookup (overload-resolution) is done at run time. Later, we can introduce a compile model to effect lookup at compile time. *)


(* ::Section:: *)
(*IEnumerator and IEnumerable*)


(* ::Input:: *)
(*ClearAll[*)
(*IEnumerable,*)
(*IEnumerator,*)
(*IEnumerableType,*)
(*IEnumeratorType,*)
(*GetEnumerator,*)
(*MoveNext,*)
(*Current];*)


(* ::Subsection:: *)
(*Abstract Types*)


(* ::Text:: *)
(*Here are the abstract specifications of the types IEnumerable and IEnumerator. This is not a full type system, but it is a start in the right direction. You can check these types by calling ProvidesTypeQ[obj,type]. *)


(* ::Input:: *)
(*IEnumerableType={GetEnumerator[]};*)


(* ::Input:: *)
(*IEnumeratorType={MoveNext[],Current};*)


(* ::Subsection:: *)
(*Contracts*)


(* ::Text:: *)
(*Unlike C#, where interfaces are abstract, our IEnumerator provide an implementation that enforces its contract by requiring every provider object -- every object that implements the interface -- to implement a private protocol, represented by MoveNext[_Integer] and Current[_Integer]. This is better than simply hoping that each provider implements the semantics properly, and it's similar to the discipline required in Scala. The private contract for the private protocol is documented in the comments of the code below.*)


(* ::Input:: *)
(*privateProtocol={MoveNext[_Integer],Current[_Integer]};*)


(* ::Text:: *)
(*To get the IEnumerator interface for an object named this, "call" the following "function" on an object (i.e., dispatch list of rewrite rules) to get a new object that implements (provides rules for) the enumerator's state machine. The new object will call the private protocol of the old object by imposing the old object's rewrite rules on expressions of the form MoveNext[i_Integer] and Current[i_Integer] that must be provided by the old object.*)


(* ::Item:: *)
(*SIDEBAR: The quotes on "call" and "function" are to remind you that what you're really doing is rewriting the expression IEnumerator[this_], after substituting an actiual object for the variable this_, into the right-hand-side of the := assignment symbol. The rewrite rule for IEnumerator is permanently installed into the global symbol table because IEnumerator is a permanent part of our programming environment. It's important to keep in mind that Mathematica implements function-calling with expression-rewriting, as can be observed by evaluating the following two expressions.  Trace[x+x/.{x->42}] rewrites first and substitutes second, i.e., rewrites x+x into 2x first, then substitutes 42 for x. Trace[Function[x,x+x][42]] substitues first and rewrites second, i.e., substitutes 42 for x and x+x for Function[x,x+x][42] in one step, and then rewrites 42+42 into 84 after the substitutions. *)


(* ::Input:: *)
(*(* This checks that its argument provides the privateProtocol. *)*)
(*IEnumerator[this_?(ProvidesTypeQ[#,privateProtocol]&)]:=*)
(*Module[{(* Variables for the state-machine. *)*)
(*i=0,*)
(*iPlus=Undefined},{*)
(*(* MoveNext is, syntactically, a method. *)*)
(*MoveNext[]:>( *)
(*(* Access this's private implementation of MoveNext[i] *)*)
(*iPlus=MoveNext[i]/.this;*)
(*(* The private protocol returns integer indices only when MoveNext stays 'inside' the sequence. *)*)
(*With[{result=(Integer===Head[iPlus])},*)
(*If[result,i=iPlus];*)
(*(* The public protocol for IEnumerator's MoveNext produces a Boolean. *)*)
(*result]),*)
(*(* Current is, syntactically, a property. *)*)
(*Current:>*)
(*(* The following effects the public, documented protocol in terms of a private protocol member, namely Current[i]. *)*)
(*If[Not[Integer===Head[iPlus]],*)
(*Throw["InvalidOperationException"],*)
(*If[i===0,*)
(*Undefined,*)
(*(* Access this's private implementation of Current[i] *)*)
(*Current[i]/.this]]*)
(*}]//Dispatch*)


(* ::Text:: *)
(*To get the IEnumerable interface for an object, call the following function, which returns a new object implementing the interface, i.e., providing rewrite rules for the interface's members.*)


(* ::Input:: *)
(*IEnumerable[this_?ObjectQ]:={*)
(*GetEnumerator[]:>IEnumerator[this]*)
(*}//Dispatch*)


(* ::Section:: *)
(*List: a Provider of IEnumerable*)


(* ::Text:: *)
(*list is a provider of IEnumerable. Its type is list; its constructor is list`list, using the class list as a namespace. It must implement the private protocol MoveNext[_Integer] and Current[_Integer]. *)


(* ::Input:: *)
(*ClearAll[list`list];*)
(*list`list[data_List]:=*)
(*Module[{len=Length[data]},*)
(*IEnumerable[{(* what follows is the required private protocol. *)*)
(*MoveNext[i_Integer]:>*)
(*With[{iPlus=i+1},*)
(*(* Required to produce an integer iff the new index is in-range. *)*)
(*If[iPlus>0&&iPlus<=len,iPlus,False]],*)
(*Current[i_Integer]:>*)
(*If[i>0&&i<=len,*)
(*data[[i]],*)
(*Throw["IndexOutOfRangeException"]]*)
(*}//Dispatch*)
(*]]*)


(* ::Section:: *)
(*ForEach*)


(* ::Text:: *)
(*Here is a straightforward implementation of forEach, which applies a function to every element of an IEnumerable for side-effect. It's a prototype for the entire suite of LINQ-ish Standard Query Operators (http://bit.ly/IS29t3). *)


(* ::Input:: *)
(*ClearAll[forEach];*)
(*forEach[enumerable_?(ProvidesTypeQ[#,IEnumerableType]&),someFunction_]:=*)
(*With[{enumerator=GetEnumerator[]/.enumerable},*)
(*While[*)
(*MoveNext[]/.enumerator,*)
(*someFunction[Current/.enumerator]*)
(*]]*)


(* ::Input:: *)
(*forEach[list`list[{"John Smith","Jim Johnson","Sue Rabon"}],Print]*)


(* ::Subsection:: *)
(*A Syntactic Improvement*)


(* ::Text:: *)
(*Leo Bushkin and I figured out how to overload Dot -- normally vector inner product -- so we can use more natural OOP notation.*)


(* ::Input:: *)
(*ClearAll[Flip];*)
(*Flip[fn_]:=Function[{x,y},fn[y,x]];*)
(*Unprotect[Dot];*)
(*SetAttributes[Dot,HoldRest];*)
(*Dot[rules_,member_]:=*)
(*Fold[Flip@ReplaceAll, List @@ rules,{Unevaluated[member]}];*)
(*Dot[rules_,member_,members__] :=((*Print@rules;Print@member;Print@members;*)*)
(*Fold[Flip@ReplaceAll, List @@ rules,{Unevaluated[member,members]}]);*)
(*Protect[Dot];*)


(* ::Text:: *)
(*Here is the more natural notation:*)


(* ::Input:: *)
(*ClearAll[forEach];*)
(*forEach[enumerable_?(ProvidesTypeQ[#,IEnumerableType]&),someFunction_]:=*)
(*With[{enumerator=enumerable.GetEnumerator[]},*)
(*While[*)
(*enumerator.MoveNext[],*)
(*someFunction[enumerator.Current]*)
(*]]*)


(* ::Input:: *)
(*forEach[list`list[{"John Smith","Jim Johnson","Sue Rabon"}],Print]*)


(* ::Text:: *)
(*Should behave well for empty lists.*)


(* ::Input:: *)
(*forEach[list`list[{}],Print]*)


(* ::Text:: *)
(*Should not reduce for things that don't provide the IEnumerable type. *)


(* ::Input:: *)
(*forEach[foobar,Print]//FullForm*)


(* ::Section:: *)
(*HashMap: a Provider of IEnumerable*)


(* ::Input:: *)
(*aHashMap=RandomSample@MapThread[Rule,{CharacterRange["a","z"],Range[26]}]//Dispatch*)


(* ::Text:: *)
(*My impl of hashMap cheats by just using the impl of list`list.*)


(* ::Input:: *)
(*hashMap`hashMap[kvs_Dispatch]:=list`list[List@@#&/@kvs[[1]]//Sort]*)


(* ::Input:: *)
(*forEach[hashMap`hashMap[aHashMap],Print]*)
