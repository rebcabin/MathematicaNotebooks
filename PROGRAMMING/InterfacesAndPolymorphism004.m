(* ::Package:: *)

(* ::Title:: *)
(*Polymorphism and Interfaces in Mathematica*)


(* ::Subtitle:: *)
(*Brian Beckman*)
(*Version of 26 December 2013*)


(* ::Section:: *)
(*Introduction*)


(* ::Text:: *)
(*Imagine the classic C# Enumerable API, leaving out Reset, and implement it for Mathematica Lists and Mathematica hash tables, which are lists of replacement rules subject to Dispatch. This serves as an illustration of a way to implement polymorphic interfaces, that is, interfaces implemented differently by different concrete providers.*)


(* ::Text:: *)
(*Then step via the now-classic duality argument [1] to interfaces for Observable and Observer.*)


(* ::Section:: *)
(*References*)


(* ::ItemNumbered:: *)
(*Enumerable Dual to Observable: http://stanford.io/1kw535m*)


(* ::ItemNumbered:: *)
(*IEnumerator Interface: http://bit.ly/1jxWqHe*)


(* ::ItemNumbered:: *)
(*Pattern-Matching Bug: http://bit.ly/1bosOYv*)


(* ::ItemNumbered:: *)
(*Standard Query Operators: http://bit.ly/IS29t3*)


(* ::ItemNumbered:: *)
(*Possible Notebook Bug: http://bit.ly/1gUoJQD*)


(* ::ItemNumbered:: *)
(*Clojure's hash-map API: http://clojure.org/data_structures*)


(* ::Section:: *)
(*Object-Oriented Programming (OOP)*)


(* ::Text:: *)
(*Represent oop's objects as explicit lists of rules that transform patterns into expressions. When an object goes out of scope, its list of rules is garbage-collected. This representation avoids stocking Mathematica's global symbol table with rules, avoiding in-the-offing explicit memory management for those rules. *)


(* ::Text:: *)
(*Every object-as-instance-of-a-class must have a rule for every class member. Likewise, an object-as-implementor-of-interfaces must have a rule for every member of every interface that the object implements. *)


(* ::Text:: *)
(*We say that an object provides implementations for its class's and interface's members, or that the object is a provider of its class's and interface's members. In both cases, the object provides implementations by providing rewrite rules. *)


(* ::Subsection:: *)
(*Type-Checking*)


(* ::Text:: *)
(*We insist that an object must be a list of rewrite rules explicitly constructed in a Dispatch table, i.e., an expression with Head Dispatch. *)


(* ::Text:: *)
(*Since Mathematica rewrites dispatches back to lists at its own discretion, we cannot be sure that every object is a Dispatch even if we insist that the user explicitly construct objects by invoking Dispatch. Instead, we can only check that an object is either a List or a Dispatch, and that's what the helper function ObjectQ does. This helper is invoked as part of the pattern-matching at call sites of functions that require objects. This invocation at the call site implements run-time type checking.*)


ClearAll[RuleQ,PatternQ,ObjectQ];

RuleQ[this_]:=With[{k=Head@this},k===Rule||k===RuleDelayed];

PatternQ[pattern_]:=True;
(* Is this even well founded question, since any expression may appear in a pattern position *)

ObjectQ[this_]:=
Switch[Head@this,
Dispatch,True,
List,(* Ensure every element is a rule *)this===Select[this,RuleQ],
_,False]


(* ::Text:: *)
(*To type-check that an object provides a certain class or interface, just check that its list-of-rules or its Dispatch contains rules for each required member. *)


(* ::Text:: *)
(*Represent types -- i.e., either classes or interfaces -- as lists of the patterns -- i.e., left-hand-sides -- required for each rule -- i.e., member -- of some concrete representation -- i.e., object. For example, every provider of IEnumerator must provide a rule with a pattern MoveNext[] -- a void-to-Boolean method -- and a pattern for Current -- a property producing a value. *)


(* ::Text:: *)
(*To check that an object provides such rules, extract the patterns on the left-hand sides of the rules of the object, then check that they are a superset of  the required patterns. A subtlety is that names must be stripped from parameters of patterns-with-parameters, leaving only the types of parameters. For example, we must convert MoveNext[i_Integer] to MoveNext[_Integer]. Accomplish this by rewriting Pattern expressions to Blank expressions (though see this known bug in the patterm-matcher\[NonBreakingSpace][3]).*)


(* ::Text:: *)
(*TODO: recurse ProvidesTypeQ on parameter types.*)


ClearAll[ProvidesTypeQ,SubsetQ,GetRules,GetPatterns,GetReplacements,StripName];

SubsetQ[A_List,B_List]:=Complement[A,B]==={};
SubsetQ[else___]:=Throw["IllegalArgumentsException: "<>ToString@{else}]

StripName[
Verbatim[Pattern][nym_,typeSpec:Verbatim[Blank][type___]]]:=typeSpec;
StripName[hd_[args___]]:=hd@@StripName/@{args};
StripName[else_]:=else

GetRules[this_?ObjectQ]:=
Sort@With[{h=Head@this},
If[h===Dispatch,
this[[1]],
If[h===List,
this,
Throw["InvalidOperationException"]]]];

GetPatterns[this_?ObjectQ]:=
StripName/@First/@GetRules@this;

GetReplacements[this_?ObjectQ]:=
#[[2]]&/@GetRules@this

ProvidesTypeQ[this_?ObjectQ,type_List]:=SubsetQ[StripName/@type,GetPatterns[this]];
ProvidesTypeQ[else___]:=Throw["IllegalArgumentsException: "<>ToString@{else}]


(* ::Text:: *)
(*This does not handle generics, type-wildcards, subtyping, and co- and contra-variance. We leave those developments for another time and place. *)


(* ::Subsubsection::Closed:: *)
(*Unit Tests*)


ProvidesTypeQ[{f[x_Integer,y_Real]:>x+y},{f[_Integer,_Real]}]


ProvidesTypeQ[{f[x_Integer,y_Real]:>x+y},{f[p_Integer,q_Real]}]


ProvidesTypeQ[{f[x_Integer,y_Real]:>x+y}//Dispatch,{f[_Integer,_Real]}]


!ProvidesTypeQ[{f[x_Integer,y_Real]:>x+y},{f[_Integer,_Integer]}]


ProvidesTypeQ[
{f[x_Integer,y_Real]:>x+y,g[s_String]:>s},{f[_Integer,_Real]}]


!ProvidesTypeQ[
{f[x_Integer,y_Real]:>x+y,g[s_String]:>s},{f[_Integer,_Real],h[_String]}]


$obj={
h[s2_String]:>s2,
f[x_Integer,y_Real]:>x+y,
g[s_String]:>s}//Sort//Dispatch;


ProvidesTypeQ[$obj,
{f[_Integer,_Real],h[_String]}]


ProvidesTypeQ[$obj,
{f[q_Integer,p_Real],h[r_String]}]


ProvidesTypeQ[$obj,{f[_Integer,_Real],h[accidentalName_String]}]


(* ::Subsection:: *)
(*Utility Functions on Types*)


(* ::Text:: *)
(*We need a way to destructively add a rule to a given object, and it's useful to have a method to check for existence of a single rule pattern; ProvidesTypeQ above checks for the existence of multiple rule patterns.*)


ClearAll[HasPattern,IsPattern,WriteRule,RemoveRule];

HasPattern[this_?ObjectQ,pattern_?RuleQ]:=
SubsetQ[{StripName@First@pattern},GetPatterns@this];
HasPattern[this_?ObjectQ,pattern_?PatternQ]:=
SubsetQ[{StripName@pattern},GetPatterns@this];

IsPattern[this_,that_?RuleQ]:=StripName@First@this===StripName@First@that;
IsPattern[this_,that_?PatternQ]:=StripName@First@this===StripName@that;

RemoveRule[this_?ObjectQ,ruleOrPattern_]:=
Sort@Select[GetRules@this,!IsPattern[#,ruleOrPattern]&];

WriteRule[this_?ObjectQ,rule_?RuleQ]:=
Append[RemoveRule[this,rule],rule]//Sort//Dispatch;
WriteRule[else___]:=Throw["IllegalArgumentsException: "<>ToString@{else}]


(* ::Subsubsection:: *)
(*Digression: Consider Clojure hash-maps*)


(* ::Text:: *)
(*Clojure has an exceptionally well conceived hash-map data structure and collection of functions. We may, in a future version of this document, emulate that structure and its API. Quoting from [6]*)


(* ::Text:: *)
(*Maps (IPersistentMap)*)
(*A Map is a collection that maps keys to values.Two different map types are provided - hashed and sorted.Hash maps require keys that correctly support hashCode and equals.Sorted maps require keys that implement Comparable, or an instance of Comparator.Hash maps provide faster access (log32N hops) vs (logN hops), but sorted maps are, well, sorted.count is O (1).conj expects another (possibly single entry) map as the item, and returns a new map which is the old map plus the entries from the new, which may overwrite entries of the old.conj also accepts a MapEntry or a vector of two items (key and value).seq returns a sequence of map entries, which are key/value pairs.Sorted map also supports rseq, which returns the entries in reverse order.Maps implement IFn, for invoke () of one argument (a key) with an optional second argument (a default value), i.e.maps are functions of their keys.nil keys and values are ok.Related functions*)
(*Create a new map : hash - map sorted - map sorted - map - by*)
(*' change' a map : assoc dissoc select - keys merge merge - with zipmap*)
(*Examine a map : get contains?find keys vals map?Examine a map entry : key val*)


(* ::Subsubsection::Closed:: *)
(*UnitTests*)


{id->Unique[]}


{id->Unique[]}.id


HasPattern[{id->Unique[]},jd]


HasPattern[{id->Unique[]},id]


HasPattern[$obj,
h[accidentalName_String]]


IsPattern[#,f[x_Integer,y_Real]]&/@GetRules@$obj==={True,False,False}


RemoveRule[$obj,
f[x_Integer,y_Real]
]===Sort@{
h[s2_String]:>s2,
g[s_String]:>s}


RemoveRule[$obj,
f[x_Integer,y_Real]:>x+y
]===Sort@{
h[s2_String]:>s2,
g[s_String]:>s}


RemoveRule[$obj,
f[z_Integer,q_Real]:>x+y
]===Sort@{
h[s2_String]:>s2,
g[s_String]:>s}


RemoveRule[$obj,
f[z_Integer,q_Real]:>z+q
]===Sort@{
h[s2_String]:>s2,
g[s_String]:>s}


RemoveRule[$obj,
j[l_List]:>Length@l]===$obj


WriteRule[$obj,j[l_List]:>Length@l]===Append[$obj,
j[l_List]:>Length@l]


WriteRule[$obj,f[q_Integer,z_Real]:>q^2+z^2]


WriteRule[$obj,f[q_Integer,z_Real]:>q^2+z^2]===Sort@{h[s2_String]:>s2,g[s_String]:>s,f[q_Integer,z_Real]:>q^2+z^2}


(* ::Subsection::Closed:: *)
(*Inheritance and Its Style of Polymorphism*)


(* ::Text:: *)
(*Particular rules could simulate inheritance and its style of polymorphism either by chasing prototype chains (upward inheritance) or dispatching at run time (downward inheritance). We don't need to take a cosmic position on this; each class and interface can do it in its own way. *)


(* ::Subsection::Closed:: *)
(*Using Objects*)


(* ::Text:: *)
(*To produce application-level results, apply objects' rules to other, arbitrarily rich work-expressions. This is "rewrite at run time." Contrast with the more familiar "rewrite at compile time," which is how standard oop systems transform work-expressions into actionable, application-level work-code. In a future version of this document, we might exhibit a compile-time option.*)


(* ::Subsection::Closed:: *)
(*Overload Resolution; Member Lookup*)


(* ::Text:: *)
(*In the current design, member lookup (overload-resolution) is done at run time. Later, we might include lookup at compile time. *)


(* ::Section:: *)
(*IEnumerator and IEnumerable*)


ClearAll[
IEnumerable,
IEnumerator,
IEnumerableType,
IEnumeratorType,
GetEnumerator,
MoveNext,
Current];


(* ::Subsection:: *)
(*Abstract Types*)


(* ::Text:: *)
(*Here are the abstract specifications of the types IEnumerable and IEnumerator. This is not a full type system, but it is a start in the right direction. You can check these types by calling ProvidesTypeQ[obj,type]. *)


ClearAll[IEnumerableType];
IEnumerableType={GetEnumerator[]};


ClearAll[IEnumeratorType];
IEnumeratorType={MoveNext[],Current};


(* ::Subsection:: *)
(*Contracts*)


(* ::Text:: *)
(*Unlike C#, where interfaces are abstract, our IEnumerator provide an implementation that enforces its contract by requiring every provider object to implement the private protocol of type MoveNext[_Integer] and Current[_Integer]. Enforcing is better than simply hoping. The private contract for the private protocol is documented in the comments of the code below.*)


ClearAll[privateProtocolType];
privateProtocolType={MoveNext[_Integer],Current[_Integer]};


(* ::Text:: *)
(*To get the IEnumerator interface for an object named this, "call" the following "function" on an object (i.e., dispatch list of rewrite rules) to get a new object that implements (provides rules for) the enumerator's state machine. The new object will call the private protocol of the old object by imposing the old object's rewrite rules on expressions of the form MoveNext[i_Integer] and Current[i_Integer] that must be provided by the old object.*)


(* ::Item:: *)
(*SIDEBAR: The quotes on "call" and "function" are to remind you that what you're really doing is rewriting the expression IEnumerator[this_], after substituting an actiual object for the variable this_, into the right-hand-side of the := assignment symbol. The rewrite rule for IEnumerator is permanently installed into the global symbol table because IEnumerator is a permanent part of our programming environment. Mathematica implements function-calling with expression-rewriting, as can be observed by evaluating the following two expressions.  Trace[x+x/.{x->42}], which rewrites first and substitutes second, i.e., rewrites x+x into 2x first, then substitutes 42 for x;  and Trace[Function[x,x+x][42]], which substitues first and rewrites second, i.e., substitutes 42 for x and x+x for Function[x,x+x][42] in one step, and then rewrites 42+42 into 84 after the substitutions. *)


ClearAll[IEnumerator];
(* Checks that the argument provides the privateProtocol. *)
IEnumerator[this_?(ProvidesTypeQ[#,privateProtocolType]&)]:=
Module[{(* Variables for the state-machine. *)
i=0,
iPlus=Undefined},{
(* MoveNext is, syntactically, a method. *)
MoveNext[]:>( 
(* Access the private implementation of MoveNext[i] *)
iPlus=MoveNext[i]/.this;
(* The private protocol returns integer indices only when MoveNext stays 'inside' the sequence. *)
With[{result=(Integer===Head[iPlus])},
If[result,i=iPlus];
(* The public protocol for IEnumerator's MoveNext produces a Boolean. *)
result]),
(* Current is, syntactically, a property. *)
Current:>
(* The following effects the public, documented protocol in terms of a private protocol member, namely Current[i]. *)
If[Not[Integer===Head[iPlus]],
Throw["InvalidOperationException"],
If[i===0,
Undefined,
(* Access this's private implementation of Current[i] *)
Current[i]/.this]]
}]//Sort//Dispatch


(* ::Text:: *)
(*To get the IEnumerable interface for an object, call the following function, which returns a new object implementing the interface, i.e., providing rewrite rules for the interface's members.*)


ClearAll[IEnumerable];
IEnumerable[this_?ObjectQ]:={
GetEnumerator[]:>IEnumerator[this]
}//Dispatch


(* ::Section:: *)
(*List: a Provider of IEnumerable*)


(* ::Text:: *)
(*list is a provider of IEnumerable. Its type is list; its constructor is list`list, using the class list as a namespace. It must implement the private protocol MoveNext[_Integer] and Current[_Integer]. *)


ClearAll[list`list];
list`list[data_List]:=
Module[{len=Length[data]},
IEnumerable[{(* what follows is the required private protocol. *)
MoveNext[i_Integer]:>
With[{iPlus=i+1},
(* Required to produce an integer iff the new index is in-range. *)
If[iPlus>0&&iPlus<=len,iPlus,False]],
Current[i_Integer]:>
If[i>0&&i<=len,
data[[i]],
Throw["IndexOutOfRangeException"]]
}//Sort//Dispatch
]]


(* ::Section:: *)
(*ForEach*)


(* ::Text:: *)
(*Here is a straightforward implementation of forEach, which applies a function to every element of an IEnumerable for side-effect. It corresponds to Mathematica's Scan. It's a prototype for the entire suite of LINQ-ish Standard Query Operators [4]. *)


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


ClearAll[Flip];
Flip[fn_]:=Function[{x,y},fn[y,x]];
Unprotect[Dot];
SetAttributes[Dot,HoldRest];
Dot[rules_,member_]:=
Fold[Flip@ReplaceAll, List @@ rules,{Unevaluated[member]}];
Dot[rules_,member_,members__] :=((*Print@rules;Print@member;Print@members;*)
Fold[Flip@ReplaceAll, List @@ rules,Unevaluated@{member,members}]);
Protect[Dot];


(* ::Subsubsection:: *)
(*Unit Test*)


(* ::Input:: *)
(*{m->{n->{p->42}}}.m.n.p*)


(* ::Text:: *)
(*Here is the more natural notation in action:*)


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
(*Should throw for things that don't provide the IEnumerable type. *)


(* ::Input:: *)
(*Catch[forEach[foobar,Print]]*)


(* ::Section:: *)
(*HashMap: a Provider of IEnumerable*)


(* ::Input:: *)
(*aHashMap=RandomSample@MapThread[Rule,{CharacterRange["a","z"],Range[26]}]//Dispatch*)


(* ::Text:: *)
(*My impl of hashMap cheats by just using the impl of list`list.*)


hashMap`hashMap[kvs_Dispatch]:=list`list[List@@#&/@kvs[[1]]//Sort]


(* ::Input:: *)
(*forEach[hashMap`hashMap[aHashMap],Print]*)


(* ::Section:: *)
(*IObserver and IObservable*)


ClearAll[
IObservable,
ISubject,
IObserver,
IObservableType,
IObserverType,
Subscribe,
OnNext, 
OnError,
OnCompleted,
Current];


(* ::Subsection:: *)
(*Delay and Force*)


(* ::Text:: *)
(*Due to a possible bug in the Notebook interface [5], we must abstract definitions of Delay and Force. It's natural to express Delay by wrapping a delayed expression in a thunk -- function of no arguments -- and then Force evaluation by calling the thunk. This is so natural that it hardly merits calling out, and I would not have done but for the fact that it sometimes does not work in the Mathematica Notebook interface, though it does work in the command-line version of Mathematica. The solution is to abstract the operations into explicit calls of Delay and Force and use Hold and ReleaseHold, which work in the Notebook and in the command-line version of Mathematica. *)


ClearAll[Delay,Force];
Delay=Hold;
Force=ReleaseHold;


(* ::Subsection:: *)
(*Abstract Types*)


ClearAll[IObservableType];
IObservableType={Subscribe[_]};


(* ::Text:: *)
(*TODO: Should be Subscribe[_IObserverType]. Also should have a type for exceptions.*)


ClearAll[IObserverType];
IObserverType={
OnNext[observation_],
OnError[exception_],
OnCompleted[]};


(* ::Subsection:: *)
(*Contracts*)


ClearAll[ISubject];
ISubject[this_?ObjectQ]:=
Module[{subscriptions={}},
{DebugReport[]:>subscriptions,
OnNext[obn_]:>
Scan[#.OnNext[obn]&,GetReplacements@subscriptions],
OnError[exc_]:>(Scan[#.OnError[exc]&,GetReplacements@subscriptions];subscriptions=Null),
OnCompleted[]:>(Scan[#.OnCompleted[]&,GetReplacements@subscriptions];subscriptions=Null),
Subscribe[that_?(ProvidesTypeQ[#,IObserverType]&)]:>
Module[{id,subscription},
id=If[HasPattern[that,SubscriptionId],
that.SubscriptionId,
Unique[]];
subscription=(id->that);
AppendTo[subscriptions,subscription];
With[{unsubscribe=Delay[subscriptions=RemoveRule[subscriptions,id]]},
unsubscribe]]
}]//Sort//Dispatch


(* ::Subsubsection:: *)
(*User-Supplied Subscription Ids*)


(* ::Input:: *)
(*ClearAll[$obr];*)
(*$obr[subscriptionId_]:={*)
(*SubscriptionId:>subscriptionId,*)
(*OnNext[msg_]:>Print[ToString@subscriptionId<>": Observer OnNext: "<> ToString@msg],*)
(*OnError[exc_]:>Print[ToString@subscriptionId<>": Observer OnError: "<> ToString@exc],*)
(*OnCompleted[]:>Print[ToString@subscriptionId<>": Observer OnCompleted!"]}//Sort//Dispatch*)


(* ::Input:: *)
(*unsubscribeFirstObserver=$myObl.Subscribe[$obr[Unique[]]]*)


(* ::Input:: *)
(*$myObl.DebugReport[]*)


(* ::Input:: *)
(*$myObl.OnNext[42]*)


(* ::Input:: *)
(*unsubscribeSecondObserver=$myObl.Subscribe[$obr[Unique[]]]*)


(* ::Input:: *)
(*$myObl.DebugReport[]*)


(* ::Input:: *)
(*$myObl.OnNext[42]*)


(* ::Input:: *)
(*Force[unsubscribeFirstObserver];*)


(* ::Input:: *)
(*$myObl.OnNext[42]*)


(* ::Subsubsection:: *)
(*System-Supplied Subscription Ids*)


(* ::Input:: *)
(*$obr[]:={*)
(*OnNext[msg_]:>Print["Observer OnNext: "<> ToString@msg],*)
(*OnError[exc_]:>Print["Observer OnError: "<> ToString@exc],*)
(*OnCompleted[]:>Print["Observer OnCompleted!"]}//Sort//Dispatch*)


(* ::Input:: *)
(*unsubscribeThirdObserver=$myObl.Subscribe[$obr[]]*)


(* ::Input:: *)
(*$myObl.DebugReport[]*)


(* ::Input:: *)
(*unsubscribeFourthObserver=$myObl.Subscribe[$obr[]]*)


(* ::Input:: *)
(*$myObl.OnNext[42]*)


(* ::Input:: *)
(*Force[unsubscribeThirdObserver]*)


(* ::Input:: *)
(*$myObl.OnNext[42]*)


(* ::Input:: *)
(*$myObl.OnCompleted[]*)


(* ::Input:: *)
(*$myObl.OnNext[42]*)


(* ::Subsection:: *)
(*GenerateWithTime*)


(* ::Input:: *)
(*ClearAll[$foobr];*)
(*$foobr[subscriptionId_]:={*)
(*SubscriptionId:>subscriptionId,*)
(*OnNext[msg_]:>($foo=ToString@subscriptionId<>": Observer OnNext: "<> ToString@msg),*)
(*OnError[exc_]:>($foo=ToString@subscriptionId<>": Observer OnError: "<> ToString@exc),*)
(*OnCompleted[]:>($foo=ToString@subscriptionId<>": Observer OnCompleted!")}//Sort//Dispatch*)


(* ::Input:: *)
(*Dynamic[$foo]*)


ClearAll[GenerateWithTime];
GenerateWithTime[
initialState_,
condition_,(* State \[Rule] Bool *)
resultSelector_, (* State \[Rule] Result *)
timeSelector_, (* State \[Rule] Real *)
iterate_ (* State \[Rule] State *)]:=
{Subscribe[obr_?(ProvidesTypeQ[#,IObserverType]&)]:>
Module[{task,getTask,
createCleanUpTask,cleanUpTask,getCleanUpTask,
state=initialState},
(* TODO: catch exceptions and call OnError *)
getTask[]:=task;
getCleanUpTask[]:=cleanUpTask;
createCleanUpTask[]:=cleanUpTask=CreateScheduledTask[
RemoveScheduledTask@getTask[];
obr.OnCompleted[];
RemoveScheduledTask@getCleanUpTask[],
{timeSelector@state}];
task=RunScheduledTask[
obr.OnNext[resultSelector@state];
state=iterate@state;
If[condition@state,
StartScheduledTask@getTask[],
StartScheduledTask@createCleanUpTask[]],
{timeSelector@state}];
]}//Sort//Dispatch;


(* ::Input:: *)
(*$foo=-1;*)
(*RemoveScheduledTask@ScheduledTasks[];*)
(*Force[*)
(*GenerateWithTime[0,#<5&,#&,(0.25)&,#+1&].Subscribe[$foobr[Unique[]]]];*)


(* ::Input:: *)
(*ScheduledTasks[]*)
