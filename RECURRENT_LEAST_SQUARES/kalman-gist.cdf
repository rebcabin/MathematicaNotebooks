(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 10.4' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1064,         20]
NotebookDataLength[    168769,       3687]
NotebookOptionsPosition[    162897,       3513]
NotebookOutlinePosition[    167307,       3622]
CellTagsIndexPosition[    167264,       3619]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Kalman Folding", "Title",
 CellChangeTimes->{{3.6698529556858463`*^9, 3.669852964981929*^9}, 
   3.6699272818174973`*^9}],

Cell["Brian Beckman", "Author",
 CellChangeTimes->{{3.669852973717079*^9, 3.6698529760933933`*^9}}],

Cell["17 April 2016", "Date",
 CellChangeTimes->{{3.6699160402711906`*^9, 3.669916050814829*^9}}],

Cell[CellGroupData[{

Cell["Preliminaries", "Chapter",
 CellChangeTimes->{{3.669844780564543*^9, 3.669844809607377*^9}}],

Cell[BoxData[{
 RowBox[{"<<", "\"\<Notation`\>\""}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"InfixNotation", "[", 
   RowBox[{
    TemplateBox[{"\[CirclePlus]"},
     "NotationTemplateTag"], ",", "Join"}], "]"}], ";"}]}], "Input",
 CellChangeTimes->{{3.669844108078168*^9, 3.669844152555797*^9}}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"ClearAll", "[", 
   RowBox[{
   "col", ",", "row", ",", "id", ",", "con", ",", "zero", ",", "len", ",", 
    "dim", ",", "m", ",", "f1", ",", "inv", ",", "pinv", ",", "scalar", ",", 
    "groundTruth"}], "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"col", "[", "xs_List", "]"}], ":=", 
   RowBox[{"List", "/@", "xs"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"row", "[", "xs_List", "]"}], ":=", 
   RowBox[{"List", "[", "xs", "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"id", "=", "IdentityMatrix"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"con", "[", 
    RowBox[{"c_", ",", "m_", ",", "n_"}], "]"}], ":=", 
   RowBox[{"ConstantArray", "[", 
    RowBox[{"c", ",", 
     RowBox[{"{", 
      RowBox[{"m", ",", "n"}], "}"}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"con", "[", 
    RowBox[{"c_", ",", "n_"}], "]"}], ":=", 
   RowBox[{"con", "[", 
    RowBox[{"c", ",", "n", ",", "n"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"zero", "[", 
    RowBox[{"m_", ",", "n_"}], "]"}], ":=", 
   RowBox[{"con", "[", 
    RowBox[{"0", ",", "m", ",", "n"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"zero", "[", "n_", "]"}], ":=", 
   RowBox[{"con", "[", 
    RowBox[{"0", ",", "n"}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"len", "=", "Length"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"dim", "[", "squareMatrix_List", "]"}], ":=", 
   RowBox[{"len", "[", 
    RowBox[{
    "squareMatrix", "\[LeftDoubleBracket]", "1", "\[RightDoubleBracket]"}], 
    "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"m", "=", "MatrixForm"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"f1", "=", 
   RowBox[{
    RowBox[{"Flatten", "[", 
     RowBox[{"#", ",", "1"}], "]"}], "&"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"inv", "=", "Inverse"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"pinv", "=", "PseudoInverse"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"scalar", "[", "m1x1_", "]"}], ":=", 
   RowBox[{"m1x1", "\[LeftDoubleBracket]", 
    RowBox[{"1", ",", "1"}], "\[RightDoubleBracket]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"groundTruth", "=", 
   RowBox[{"(", GridBox[{
      {
       RowBox[{"-", "3"}]},
      {"9"},
      {
       RowBox[{"-", "4"}]},
      {
       RowBox[{"-", "5"}]}
     }], ")"}]}], ";"}]}], "Input",
 CellChangeTimes->{
  3.629812740599895*^9, {3.629846132870491*^9, 3.629846136780778*^9}, {
   3.631926926338415*^9, 3.631926934992784*^9}, {3.669627854288382*^9, 
   3.669627894520998*^9}, {3.669822874162187*^9, 3.6698229917431593`*^9}, {
   3.669831141825832*^9, 3.66983115083389*^9}, {3.669836766230214*^9, 
   3.669836767810321*^9}, {3.669837356988452*^9, 3.669837370962223*^9}, {
   3.66984001115753*^9, 3.66984002335413*^9}, {3.669843500765081*^9, 
   3.669843538402874*^9}, {3.669843793723061*^9, 3.669843831038666*^9}, {
   3.669844381796507*^9, 3.669844525523814*^9}, 3.6698445679860277`*^9, {
   3.669844703591049*^9, 3.6698447308692293`*^9}, {3.669850414763523*^9, 
   3.669850469016633*^9}, {3.6698505009458303`*^9, 3.669850519560555*^9}, {
   3.669853493949973*^9, 3.6698535448891478`*^9}}]
}, Open  ]],

Cell[CellGroupData[{

Cell["Time-Independent Model", "Chapter",
 CellChangeTimes->{{3.6024270310824833`*^9, 3.602427034102489*^9}, 
   3.629844953075777*^9, {3.669596552029359*^9, 3.66959656266508*^9}, {
   3.669916057134766*^9, 3.6699160652141933`*^9}}],

Cell["\<\
Unit observation covariance, Generalized below in the time-dependent model.\
\>", "Text",
 CellChangeTimes->{{3.66991607021412*^9, 3.669916079549821*^9}, {
  3.669927293126544*^9, 3.6699272999263887`*^9}}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"ClearAll", "[", "kalman", "]"}], ";"}], "\n", 
 RowBox[{
  RowBox[{
   RowBox[{"kalman", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"x_", ",", "P_"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"A_", ",", "z_"}], "}"}]}], "]"}], ":=", "\[IndentingNewLine]", 
   
   RowBox[{"Module", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"D", ",", "K"}], "}"}], ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"D", "=", 
       RowBox[{
        RowBox[{"id", "[", 
         RowBox[{"len", "[", "z", "]"}], "]"}], "+", 
        RowBox[{"A", ".", "P", ".", 
         SuperscriptBox["A", "\[Transpose]"]}]}]}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"K", "=", 
       RowBox[{"P", ".", 
        RowBox[{"A", "\[Transpose]"}], ".", 
        RowBox[{"inv", "[", "D", "]"}]}]}], ";", "\[IndentingNewLine]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"x", "+", 
         RowBox[{"K", ".", 
          RowBox[{"(", 
           RowBox[{"z", "-", 
            RowBox[{"A", ".", "x"}]}], ")"}]}]}], ",", 
        RowBox[{"P", "-", 
         RowBox[{"K", ".", "D", ".", 
          SuperscriptBox["K", "\[Transpose]"]}]}]}], "}"}]}]}], "]"}]}], 
  ";"}]}], "Input",
 CellChangeTimes->{{3.668113230855397*^9, 3.66811326296312*^9}, {
   3.668113313392239*^9, 3.66811336747112*^9}, {3.6681140352071657`*^9, 
   3.668114038413567*^9}, {3.668114705772566*^9, 3.668114794519589*^9}, 
   3.6681148270239487`*^9, {3.668114861213798*^9, 3.6681148631421747`*^9}, {
   3.668115024073289*^9, 3.6681150246304693`*^9}, {3.6681170057833357`*^9, 
   3.6681170073086443`*^9}, {3.668258334294606*^9, 3.668258351947137*^9}, {
   3.6696756441589537`*^9, 3.669675679055326*^9}, {3.669822892707678*^9, 
   3.669822893257841*^9}, {3.6698229697444563`*^9, 3.669822974462617*^9}, {
   3.66983340716638*^9, 3.669833407373623*^9}, 3.669840039210088*^9, 
   3.66984165619131*^9, 3.6699126687159367`*^9, {3.669915973156136*^9, 
   3.669915973832649*^9}}],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"ClearAll", "[", "testCase", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"m", "/@", 
  RowBox[{"(", 
   RowBox[{"testCase", "=", 
    RowBox[{"{", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"{", 
          RowBox[{"1", ",", "0.", ",", "0.", ",", "0."}], "}"}], "}"}], ",", 
        RowBox[{"{", 
         RowBox[{"-", "2.28442"}], "}"}]}], "}"}], ",", "\[IndentingNewLine]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"{", 
          RowBox[{"1", ",", "1.", ",", "1.", ",", "1."}], "}"}], "}"}], ",", 
        RowBox[{"{", 
         RowBox[{"-", "4.83168"}], "}"}]}], "}"}], ",", "\[IndentingNewLine]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"{", 
          RowBox[{"1", ",", 
           RowBox[{"-", "1."}], ",", "1.", ",", 
           RowBox[{"-", "1."}]}], "}"}], "}"}], ",", 
        RowBox[{"{", 
         RowBox[{"-", "10.4601"}], "}"}]}], "}"}], ",", "\[IndentingNewLine]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"{", 
          RowBox[{"1", ",", 
           RowBox[{"-", "2."}], ",", "4.", ",", 
           RowBox[{"-", "8."}]}], "}"}], "}"}], ",", 
        RowBox[{"{", "1.40488", "}"}]}], "}"}], ",", "\[IndentingNewLine]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"{", 
          RowBox[{"1", ",", "2.", ",", "4.", ",", "8."}], "}"}], "}"}], ",", 
        RowBox[{"{", 
         RowBox[{"-", "40.8079"}], "}"}]}], "}"}]}], "}"}]}], 
   ")"}]}]}], "Input",
 CellChangeTimes->{{3.669836772187957*^9, 3.669836833247377*^9}, {
  3.669844685894539*^9, 3.6698446970945473`*^9}}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"{", 
         RowBox[{"1", ",", "0.`", ",", "0.`", ",", "0.`"}], "}"}]},
       {
        RowBox[{"-", "2.28442`"}]}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]], ",", 
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"{", 
         RowBox[{"1", ",", "1.`", ",", "1.`", ",", "1.`"}], "}"}]},
       {
        RowBox[{"-", "4.83168`"}]}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]], ",", 
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"{", 
         RowBox[{"1", ",", 
          RowBox[{"-", "1.`"}], ",", "1.`", ",", 
          RowBox[{"-", "1.`"}]}], "}"}]},
       {
        RowBox[{"-", "10.4601`"}]}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]], ",", 
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"{", 
         RowBox[{"1", ",", 
          RowBox[{"-", "2.`"}], ",", "4.`", ",", 
          RowBox[{"-", "8.`"}]}], "}"}]},
       {"1.40488`"}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]], ",", 
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"{", 
         RowBox[{"1", ",", "2.`", ",", "4.`", ",", "8.`"}], "}"}]},
       {
        RowBox[{"-", "40.8079`"}]}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]]}], "}"}]], "Output",
 CellChangeTimes->{{3.669836797162375*^9, 3.669836850945941*^9}, 
   3.669840041195919*^9, 3.669840205379318*^9, 3.669840630067752*^9, 
   3.669841664815613*^9, 3.669843550246779*^9, 3.669843849888214*^9, 
   3.669844697833145*^9, 3.6698505361652822`*^9, 3.669851964061302*^9, 
   3.669852994092271*^9, 3.66985355691544*^9, 3.669853998758136*^9, 
   3.669890989942973*^9, 3.669896284293085*^9, 3.6699176353267937`*^9, 
   3.669918502739876*^9, {3.669928977894611*^9, 3.669929007047667*^9}}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"m", "@", 
  RowBox[{"(", 
   RowBox[{
    RowBox[{
     RowBox[{"m", "/@", "#"}], "&"}], "/@", "\[IndentingNewLine]", 
    RowBox[{"Chop", "[", 
     RowBox[{"FoldList", "[", "\[IndentingNewLine]", 
      RowBox[{"kalman", ",", "\[IndentingNewLine]", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"col", "[", 
          RowBox[{"{", 
           RowBox[{"0", ",", "0", ",", "0", ",", "0"}], "}"}], "]"}], ",", 
         "\[IndentingNewLine]", 
         RowBox[{
          RowBox[{"id", "[", "4", "]"}], "*", "1000.0"}]}], "}"}], ",", 
       "\[IndentingNewLine]", "testCase"}], "\[IndentingNewLine]", "]"}], 
     "]"}]}], ")"}]}]], "Input",
 CellChangeTimes->{{3.6695963875782537`*^9, 3.669596388137123*^9}, 
   3.669627906374621*^9, {3.6696755514347277`*^9, 3.669675553936233*^9}, {
   3.669823001511177*^9, 3.669823006973716*^9}, {3.669832934096292*^9, 
   3.6698330177812223`*^9}, {3.669833201457459*^9, 3.66983320302964*^9}, {
   3.669833258756941*^9, 3.669833299417728*^9}, 3.669836782213491*^9, {
   3.6698368395348263`*^9, 3.669836841133628*^9}, {3.669851952873661*^9, 
   3.669851953238277*^9}, {3.669851990556367*^9, 3.669852011163945*^9}, 
   3.669890928153307*^9, 3.669915975834077*^9}],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {"0"},
          {"0"},
          {"0"},
          {"0"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]], 
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {"1000.`", "0", "0", "0"},
          {"0", "1000.`", "0", "0"},
          {"0", "0", "1000.`", "0"},
          {"0", "0", "0", "1000.`"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]]},
     {
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {
           RowBox[{"-", "2.282137862137862`"}]},
          {"0"},
          {"0"},
          {"0"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]], 
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {"0.9990009990009412`", "0", "0", "0"},
          {"0", "1000.`", "0", "0"},
          {"0", "0", "1000.`", "0"},
          {"0", "0", "0", "1000.`"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]]},
     {
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {
           RowBox[{"-", "2.282986295179269`"}]},
          {
           RowBox[{"-", "0.8492814744487608`"}]},
          {
           RowBox[{"-", "0.8492814744487608`"}]},
          {
           RowBox[{"-", "0.8492814744487608`"}]}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]], 
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {"0.998668552855665`", 
           RowBox[{"-", "0.3327785914214145`"}], 
           RowBox[{"-", "0.3327785914214145`"}], 
           RowBox[{"-", "0.3327785914214145`"}]},
          {
           RowBox[{"-", "0.3327785914214145`"}], "666.8886299871448`", 
           RowBox[{"-", "333.1113700128552`"}], 
           RowBox[{"-", "333.1113700128552`"}]},
          {
           RowBox[{"-", "0.3327785914214145`"}], 
           RowBox[{"-", "333.1113700128552`"}], "666.8886299871448`", 
           RowBox[{"-", "333.1113700128552`"}]},
          {
           RowBox[{"-", "0.3327785914214145`"}], 
           RowBox[{"-", "333.1113700128552`"}], 
           RowBox[{"-", "333.1113700128552`"}], "666.8886299871448`"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]]},
     {
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {
           RowBox[{"-", "2.2874882356667148`"}]},
          {"1.406753311672082`"},
          {
           RowBox[{"-", "5.355723902382095`"}]},
          {"1.406753311672082`"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]], 
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {"0.9980044897732641`", "0", 
           RowBox[{"-", "0.9975057369048116`"}], "0"},
          {"0", "500.12496875781056`", "0", 
           RowBox[{"-", "499.87503124218944`"}]},
          {
           RowBox[{"-", "0.9975057369048115`"}], "0", "1.4967573582258638`", 
           "0"},
          {"0", 
           RowBox[{"-", "499.87503124218944`"}], "0", "500.12496875781056`"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]]},
     {
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {
           RowBox[{"-", "2.2939903233863417`"}]},
          {"7.923470322389598`"},
          {
           RowBox[{"-", "5.3448809476317765`"}]},
          {
           RowBox[{"-", "5.1153952018048265`"}]}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]], 
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {"0.9975079870225654`", "0.49761985086406924`", 
           RowBox[{"-", "0.9966777633228829`"}], 
           RowBox[{"-", "0.49803460319851184`"}]},
          {"0.4976198508640693`", "1.3855043227356418`", 
           RowBox[{"-", "0.8298364707526127`"}], 
           RowBox[{"-", "0.7198813058349174`"}]},
          {
           RowBox[{"-", "0.9966777633228828`"}], 
           RowBox[{"-", "0.8298364707526126`"}], "1.495376620150582`", 
           "0.8305281164191844`"},
          {
           RowBox[{"-", "0.49803460319851184`"}], 
           RowBox[{"-", "0.7198813058349742`"}], "0.8305281164191844`", 
           "0.5537868578484222`"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]]},
     {
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {
           RowBox[{"-", "2.974226591552785`"}]},
          {"7.262403743565966`"},
          {
           RowBox[{"-", "4.210511281564406`"}]},
          {
           RowBox[{"-", "4.453777642335395`"}]}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]], 
      TagBox[
       RowBox[{"(", "\[NoBreak]", GridBox[{
          {"0.48545809498345005`", "0", 
           RowBox[{"-", "0.14277759330119488`"}], "0"},
          {"0", "0.9019078610942634`", "0", 
           RowBox[{"-", "0.2358817799647956`"}]},
          {
           RowBox[{"-", "0.14277759330119477`"}], "0", "0.07140307440990923`",
            "0"},
          {"0", 
           RowBox[{"-", "0.23588177996485243`"}], "0", 
           "0.06938393180670932`"}
         },
         GridBoxAlignment->{
          "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
           "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
         GridBoxSpacings->{"Columns" -> {
             Offset[0.27999999999999997`], {
              Offset[0.7]}, 
             Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
             Offset[0.2], {
              Offset[0.4]}, 
             Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
       Function[BoxForm`e$, 
        MatrixForm[BoxForm`e$]]]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output",
 CellChangeTimes->{{3.669596313263074*^9, 3.669596317823225*^9}, {
   3.6695963808945427`*^9, 3.6695963896905212`*^9}, 3.669627314619377*^9, 
   3.669627922729704*^9, 3.669675899429034*^9, 3.669822825597576*^9, 
   3.6698230075328836`*^9, {3.669832937136129*^9, 3.66983301852802*^9}, 
   3.6698332144142714`*^9, {3.6698332599614973`*^9, 3.6698333077144938`*^9}, 
   3.6698334103854094`*^9, {3.66983684229504*^9, 3.6698368510132523`*^9}, 
   3.669840052498838*^9, 3.669840205448661*^9, 3.669840630137002*^9, 
   3.6698416745746*^9, 3.669843592489896*^9, 3.669843852279305*^9, 
   3.6698447507683764`*^9, 3.669850538024596*^9, {3.669851954218512*^9, 
   3.669852011883193*^9}, 3.669852995289073*^9, 3.6698535569919786`*^9, 
   3.6698539987896767`*^9, 3.669890928860914*^9, 3.669890990003372*^9, 
   3.669896287025753*^9, 3.669917635395647*^9, 3.669918502808271*^9, 
   3.66992897793461*^9, 3.6699290085258017`*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell["Insensitivity to A-Priori Observation Covariance", "Subchapter",
 CellChangeTimes->{{3.669890956363908*^9, 3.6698909764650183`*^9}}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"m", "/@", 
  RowBox[{"Chop", "[", 
   RowBox[{"Fold", "[", 
    RowBox[{"kalman", ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"col", "[", 
        RowBox[{"{", 
         RowBox[{"0", ",", "0", ",", "0", ",", "0"}], "}"}], "]"}], ",", 
       RowBox[{
        RowBox[{"id", "[", "4", "]"}], "*", "1000.0"}]}], "}"}], ",", 
     "testCase"}], "]"}], "]"}]}]], "Input",
 CellChangeTimes->{{3.6695963875782537`*^9, 3.669596388137123*^9}, 
   3.669627906374621*^9, {3.6696755514347277`*^9, 3.669675553936233*^9}, {
   3.669823001511177*^9, 3.669823006973716*^9}, {3.669832934096292*^9, 
   3.6698330177812223`*^9}, {3.669833201457459*^9, 3.66983320302964*^9}, {
   3.669833258756941*^9, 3.669833299417728*^9}, 3.669836782213491*^9, {
   3.6698368395348263`*^9, 3.669836841133628*^9}, {3.669851952873661*^9, 
   3.669851953238277*^9}, {3.669851990556367*^9, 3.669852061116693*^9}, 
   3.669915978516582*^9}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"-", "2.974226591552785`"}]},
       {"7.262403743565966`"},
       {
        RowBox[{"-", "4.210511281564406`"}]},
       {
        RowBox[{"-", "4.453777642335395`"}]}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]], ",", 
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {"0.48545809498345005`", "0", 
        RowBox[{"-", "0.14277759330119488`"}], "0"},
       {"0", "0.9019078610942634`", "0", 
        RowBox[{"-", "0.2358817799647956`"}]},
       {
        RowBox[{"-", "0.14277759330119477`"}], "0", "0.07140307440990923`", 
        "0"},
       {"0", 
        RowBox[{"-", "0.23588177996485243`"}], "0", "0.06938393180670932`"}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]]}], "}"}]], "Output",
 CellChangeTimes->{{3.669852039801221*^9, 3.6698520470761433`*^9}, 
   3.669852997082272*^9, 3.669853557057708*^9, 3.669853998849139*^9, 
   3.66989629503328*^9, 3.6699176354898033`*^9, 3.6699185029074707`*^9, 
   3.669928978031938*^9, 3.669929011002204*^9}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"m", "/@", 
  RowBox[{"Chop", "[", 
   RowBox[{"Fold", "[", 
    RowBox[{"kalman", ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"col", "[", 
        RowBox[{"{", 
         RowBox[{"0", ",", "0", ",", "0", ",", "0"}], "}"}], "]"}], ",", 
       RowBox[{
        RowBox[{"id", "[", "4", "]"}], "*", "1000000.0"}]}], "}"}], ",", 
     "testCase"}], "]"}], "]"}]}]], "Input",
 CellChangeTimes->{{3.6695963875782537`*^9, 3.669596388137123*^9}, 
   3.669627906374621*^9, {3.6696755514347277`*^9, 3.669675553936233*^9}, {
   3.669823001511177*^9, 3.669823006973716*^9}, {3.669832934096292*^9, 
   3.6698330177812223`*^9}, {3.669833201457459*^9, 3.66983320302964*^9}, {
   3.669833258756941*^9, 3.669833299417728*^9}, 3.669836782213491*^9, {
   3.6698368395348263`*^9, 3.669836841133628*^9}, {3.669851952873661*^9, 
   3.669851953238277*^9}, {3.669851990556367*^9, 3.669852072889407*^9}, 
   3.669915981625524*^9}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"-", "2.975068870837714`"}]},
       {"7.2700040514448965`"},
       {
        RowBox[{"-", "4.210387267080661`"}]},
       {
        RowBox[{"-", "4.4557996407308424`"}]}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]], ",", 
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {"0.4857140292515405`", "0", 
        RowBox[{"-", "0.14285706321444525`"}], "0"},
       {"0", "0.9027769068818003`", "0", 
        RowBox[{"-", "0.23611088155041637`"}]},
       {
        RowBox[{"-", "0.14285706321444525`"}], "1.0540934791691825`*^-10", 
        "0.07142854587375425`", "0"},
       {"0", 
        RowBox[{"-", "0.23611088160862403`"}], "0", "0.06944438392173413`"}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]]}], "}"}]], "Output",
 CellChangeTimes->{3.669852073902739*^9, 3.669852998282872*^9, 
  3.6698535571248837`*^9, 3.669853998915626*^9, 3.669896295115479*^9, 
  3.669917635565507*^9, 3.669918502971589*^9, 3.669928978095278*^9, 
  3.669929012220892*^9}]
}, Open  ]]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Time-Dependent Filter", "Chapter",
 CellChangeTimes->{{3.669596919924898*^9, 3.669596924636051*^9}}],

Cell[CellGroupData[{

Cell["Process-Noise Matrix", "Subchapter",
 CellChangeTimes->{{3.669893088798563*^9, 3.669893093245285*^9}}],

Cell["\<\
Process noise always enters the last state, and with a constant standard \
deviation. We must integrate the process-noise covariance matrix over a \
discrete time interval and a particular model of system dynamics to add its \
effects to the filter.\
\>", "Text",
 CellChangeTimes->{{3.669895118788309*^9, 3.669895150355185*^9}, {
  3.669895190843395*^9, 3.669895231399893*^9}, {3.669916733276064*^9, 
  3.669916751947448*^9}}],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"ClearAll", "[", 
   RowBox[{"Q", ",", "\[CapitalPhi]", ",", "\[Sigma]x", ",", "t"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Q", "=", 
   RowBox[{"(", GridBox[{
      {"0", "0", "0"},
      {"0", "0", "0"},
      {"0", "0", "1"}
     }], ")"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"\[CapitalPhi]", "[", "t_", "]"}], ":=", 
   RowBox[{"(", GridBox[{
      {"1", "t", 
       FractionBox[
        SuperscriptBox["t", "2"], "2"]},
      {"0", "1", "t"},
      {"0", "0", "1"}
     }], ")"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Integrate", "[", 
   RowBox[{
    RowBox[{
     RowBox[{"\[CapitalPhi]", "[", "t", "]"}], ".", "Q", ".", 
     SuperscriptBox[
      RowBox[{"\[CapitalPhi]", "[", "t", "]"}], "\[Transpose]"]}], ",", 
    RowBox[{"{", 
     RowBox[{"t", ",", "0", ",", "dt"}], "}"}]}], "]"}], "//", 
  "m"}]}], "Input",
 CellChangeTimes->{{3.6698934467876453`*^9, 3.6698935759559507`*^9}, {
  3.669893665696396*^9, 3.66989367778327*^9}, {3.669893737622409*^9, 
  3.669893817344388*^9}, {3.669893877166535*^9, 3.669894024737306*^9}, {
  3.6698993446085243`*^9, 3.6698993446700163`*^9}}],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      FractionBox[
       SuperscriptBox["dt", "5"], "20"], 
      FractionBox[
       SuperscriptBox["dt", "4"], "8"], 
      FractionBox[
       SuperscriptBox["dt", "3"], "6"]},
     {
      FractionBox[
       SuperscriptBox["dt", "4"], "8"], 
      FractionBox[
       SuperscriptBox["dt", "3"], "3"], 
      FractionBox[
       SuperscriptBox["dt", "2"], "2"]},
     {
      FractionBox[
       SuperscriptBox["dt", "3"], "6"], 
      FractionBox[
       SuperscriptBox["dt", "2"], "2"], "dt"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output",
 CellChangeTimes->{
  3.6698939322999573`*^9, {3.669893977306632*^9, 3.669894029442347*^9}, 
   3.6698963042177877`*^9, 3.6698993463554993`*^9, 3.6698997049473953`*^9, 
   3.6699176368885803`*^9, 3.669918504308343*^9, 3.669928978174108*^9, 
   3.6699290148061457`*^9}]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["As A Fold", "Subchapter",
 CellChangeTimes->{{3.669893272822137*^9, 3.669893280142438*^9}}],

Cell[CellGroupData[{

Cell[TextData[{
 StyleBox["gen", "Code"],
 " is a way to generate random samples from a covariance matrix. This only \
honors the diagonal terms, but each diagonal term individually."
}], "Item",
 CellChangeTimes->{{3.669892970267221*^9, 3.66989302113589*^9}, {
  3.669895086916964*^9, 3.669895098796315*^9}, {3.669895686477333*^9, 
  3.669895687656653*^9}, {3.669896328851823*^9, 3.6698963332995462`*^9}}],

Cell[TextData[{
 StyleBox["kalman", "Code"],
 " is the Kalman filter with time-dependent dynamical process (also called ",
 StyleBox["state estimation",
  FontSlant->"Italic"],
 ")."
}], "Item",
 CellChangeTimes->{{3.669892970267221*^9, 3.669893024431725*^9}, {
   3.669893055894322*^9, 3.66989305951081*^9}, {3.669895718531981*^9, 
   3.669895723219805*^9}, {3.6698963416591578`*^9, 3.669896373737288*^9}, {
   3.669896508012918*^9, 3.669896511924958*^9}, 3.6699159929205637`*^9}],

Cell[TextData[{
 "At the bottom, we validate that ",
 StyleBox["kalman", "Code"],
 " degenerates to the original static case (also called ",
 StyleBox["parameter estimation",
  FontSlant->"Italic"],
 ")."
}], "Item",
 CellChangeTimes->{{3.669892970267221*^9, 3.669893067270576*^9}, {
   3.669895329781928*^9, 3.6698953297926493`*^9}, {3.6698957401952744`*^9, 
   3.66989575231469*^9}, {3.669896354211443*^9, 3.669896355842102*^9}, {
   3.669896521116231*^9, 3.669896522980103*^9}, 3.669915996037307*^9, {
   3.669917246223737*^9, 3.669917254503613*^9}}],

Cell[TextData[{
 Cell[BoxData[
  FormBox["\[CapitalXi]", TraditionalForm]]],
 " is the integrated process-noise matrix. It depend on time and on the size \
of the time step."
}], "Item",
 CellChangeTimes->{{3.669892970267221*^9, 3.669893067270576*^9}, {
  3.669895329781928*^9, 3.669895379058209*^9}, {3.6698955841085663`*^9, 
  3.669895615161517*^9}, {3.669895678477315*^9, 3.669895678766088*^9}, {
  3.6698957560423183`*^9, 3.669895764306785*^9}, {3.669916785921588*^9, 
  3.669916798586093*^9}}],

Cell[TextData[{
 "\[CapitalPhi] is the integral of the system-dynamics matrix F; more \
precisely, it is ",
 StyleBox["Exp[F t]", "Input"],
 StyleBox[". ", "Text"]
}], "Item",
 CellChangeTimes->{{3.669892970267221*^9, 3.669893067270576*^9}, {
  3.669895329781928*^9, 3.669895379058209*^9}, {3.6698955841085663`*^9, 
  3.669895615161517*^9}, {3.669895678477315*^9, 3.669895678766088*^9}, {
  3.6698957560423183`*^9, 3.669895764306785*^9}, {3.669916785921588*^9, 
  3.6699168220000896`*^9}, {3.6699171510687847`*^9, 3.669917180462338*^9}}],

Cell[TextData[StyleBox["\[CapitalGamma] is time-step integral of system \
response G propagated by \[CapitalPhi].", "Text"]], "Item",
 CellChangeTimes->{{3.669892970267221*^9, 3.669893067270576*^9}, {
  3.669895329781928*^9, 3.669895379058209*^9}, {3.6698955841085663`*^9, 
  3.669895615161517*^9}, {3.669895678477315*^9, 3.669895678766088*^9}, {
  3.6698957560423183`*^9, 3.669895764306785*^9}, {3.669916785921588*^9, 
  3.6699168220000896`*^9}, {3.6699171510687847`*^9, 3.66991722084899*^9}}],

Cell[TextData[{
 Cell[BoxData[
  FormBox["\[CapitalXi]", TraditionalForm]]],
 ", ",
 Cell[BoxData[
  FormBox["\[CapitalPhi]", TraditionalForm]]],
 ", ",
 Cell[BoxData[
  FormBox["\[CapitalGamma]", TraditionalForm]]],
 ", and ",
 Cell[BoxData[
  FormBox["u", TraditionalForm]]],
 " may all depend on time and on the size of the time step. "
}], "Item",
 CellChangeTimes->{{3.669892970267221*^9, 3.669893067270576*^9}, {
  3.669895329781928*^9, 3.669895379058209*^9}, {3.6698955841085663`*^9, 
  3.669895680973752*^9}, {3.669895768442128*^9, 3.669895768498581*^9}}]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"ClearAll", "[", 
    RowBox[{"gen", ",", "kalman"}], "]"}], ";"}], "\[IndentingNewLine]", 
  RowBox[{"(*", " ", 
   RowBox[{"generate", " ", "noisy", " ", "fake", " ", "observations"}], " ", 
   "*)"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{
    RowBox[{"gen", "[", "P_", "]"}], ":=", "\[IndentingNewLine]", 
    RowBox[{"Module", "[", 
     RowBox[{
      RowBox[{"{", "\[Sigma]s", "}"}], ",", "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"\[Sigma]s", "=", 
        RowBox[{"Sqrt", "/@", 
         RowBox[{"Table", "[", 
          RowBox[{
           RowBox[{"P", "[", 
            RowBox[{"[", 
             RowBox[{"i", ",", "i"}], "]"}], "]"}], ",", 
           RowBox[{"{", 
            RowBox[{"i", ",", 
             RowBox[{"dim", "[", "P", "]"}]}], "}"}]}], "]"}]}]}], ";", 
       "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{
         RowBox[{"If", "[", 
          RowBox[{
           RowBox[{
            RowBox[{"Chop", "[", "#", "]"}], "\[Equal]", "0"}], ",", "0", ",", 
           RowBox[{"(*", " ", 
            RowBox[{"RandomVariate", " ", 
             RowBox[{"can", "'"}], "t", " ", "handle", " ", "a", " ", "zero", 
             " ", "\[Sigma]"}], " ", "*)"}], "\[IndentingNewLine]", 
           RowBox[{"RandomVariate", "[", 
            RowBox[{"NormalDistribution", "[", 
             RowBox[{"0.0", ",", "#"}], "]"}], "]"}]}], "]"}], "&"}], "/@", 
        "\[Sigma]s"}]}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"this", " ", "is", " ", "a", " ", "function"}], "-", 
    RowBox[{
     RowBox[{"maker", " ", "--"}], " ", "give", " ", "it", " ", "an", " ", 
     "obs", " ", "noise", " ", "covar", " ", "and", " ", "get", " ", "a", " ",
      "foldable", " ", "function"}]}], " ", "*)"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{
    RowBox[{
     RowBox[{"kalman", "[", "\[CapitalZeta]_", "]"}], "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"x_", ",", "P_"}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{
       "\[CapitalXi]_", ",", "\[CapitalPhi]_", ",", "\[CapitalGamma]_", ",", 
        "u_", ",", "A_", ",", "z_"}], "}"}]}], "]"}], ":=", 
    "\[IndentingNewLine]", 
    RowBox[{"Module", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"x2", ",", "P2", ",", "D", ",", "K"}], "}"}], ",", 
      "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"x2", "=", 
        RowBox[{
         RowBox[{"\[CapitalPhi]", ".", "x"}], "+", 
         RowBox[{"\[CapitalGamma]", ".", "u"}]}]}], ";", 
       "\[IndentingNewLine]", 
       RowBox[{"P2", "=", 
        RowBox[{"\[CapitalXi]", "+", 
         RowBox[{"\[CapitalPhi]", ".", "P", ".", 
          SuperscriptBox["\[CapitalPhi]", "\[Transpose]"]}]}]}], ";", 
       "\[IndentingNewLine]", 
       RowBox[{"(*", " ", 
        RowBox[{
         RowBox[{
         "below", " ", "this", " ", "line", " ", "is", " ", "just", " ", 
          "like", " ", "the", " ", "time"}], "-", 
         RowBox[{"independent", " ", 
          RowBox[{"model", "!"}]}]}], " ", "*)"}], "\[IndentingNewLine]", 
       RowBox[{"D", "=", 
        RowBox[{"\[CapitalZeta]", "+", 
         RowBox[{"A", ".", "P2", ".", 
          SuperscriptBox["A", "\[Transpose]"]}]}]}], ";", 
       "\[IndentingNewLine]", 
       RowBox[{"K", "=", 
        RowBox[{"P2", ".", 
         SuperscriptBox["A", "\[Transpose]"], ".", 
         RowBox[{"inv", "[", "D", "]"}]}]}], ";", "\[IndentingNewLine]", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"x2", "+", 
          RowBox[{"K", ".", 
           RowBox[{"(", 
            RowBox[{"z", "-", 
             RowBox[{"A", ".", "x2"}]}], ")"}]}]}], ",", 
         RowBox[{"P2", "-", 
          RowBox[{"K", ".", "D", ".", 
           SuperscriptBox["K", "\[Transpose]"]}]}]}], "}"}]}]}], "]"}]}], 
   ";"}], "\n", 
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
    "show", " ", "that", " ", "it", " ", "degenerates", " ", "to", " ", "the",
      " ", "time"}], "-", 
    RowBox[{"independent", " ", "case"}]}], " ", 
   "*)"}]}], "\[IndentingNewLine]", 
 RowBox[{"m", "/@", "\[IndentingNewLine]", 
  RowBox[{"Chop", "[", "\[IndentingNewLine]", 
   RowBox[{"With", "[", 
    RowBox[{
     RowBox[{"{", "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"\[CapitalXi]", "=", 
        RowBox[{"zero", "[", "4", "]"}]}], ",", 
       RowBox[{"\[CapitalZeta]", "=", 
        RowBox[{"id", "[", "1", "]"}]}], ",", "\[IndentingNewLine]", 
       RowBox[{"\[CapitalPhi]", "=", 
        RowBox[{"id", "[", "4", "]"}]}], ",", 
       RowBox[{"\[CapitalGamma]", "=", 
        RowBox[{"zero", "[", 
         RowBox[{"4", ",", "1"}], "]"}]}], ",", 
       RowBox[{"u", "=", 
        RowBox[{"zero", "[", "1", "]"}]}]}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Fold", "[", "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"kalman", "[", "\[CapitalZeta]", "]"}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"col", "[", 
          RowBox[{"{", 
           RowBox[{"0", ",", "0", ",", "0", ",", "0"}], "}"}], "]"}], ",", 
         "\[IndentingNewLine]", 
         RowBox[{
          RowBox[{"id", "[", "4", "]"}], "*", "1000.0"}]}], "}"}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{
         RowBox[{
          RowBox[{"{", 
           RowBox[{
           "\[CapitalXi]", ",", "\[CapitalPhi]", ",", "\[CapitalGamma]", ",", 
            "u"}], "}"}], "\[CirclePlus]", "#"}], "&"}], "/@", "testCase"}]}],
       "\[IndentingNewLine]", "]"}]}], "]"}], "]"}]}]}], "Input",
 CellChangeTimes->{{3.6698530457670403`*^9, 3.669853092865712*^9}, {
   3.669853991896016*^9, 3.669854011949581*^9}, {3.669892827186557*^9, 
   3.66989284225704*^9}, 3.669895070468794*^9, 3.669895252863975*^9, {
   3.669895474567317*^9, 3.6698954838068666`*^9}, {3.669895518017056*^9, 
   3.669895540612694*^9}, 3.669896567921982*^9, {3.669896827258924*^9, 
   3.669896852824854*^9}, 3.669900017582221*^9, 3.6699000497907457`*^9, {
   3.669915998412506*^9, 3.669916002252948*^9}, {3.669927334832863*^9, 
   3.669927440786106*^9}}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {
        RowBox[{"-", "2.974226591552785`"}]},
       {"7.262403743565966`"},
       {
        RowBox[{"-", "4.210511281564406`"}]},
       {
        RowBox[{"-", "4.453777642335395`"}]}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]], ",", 
   TagBox[
    RowBox[{"(", "\[NoBreak]", GridBox[{
       {"0.48545809498345005`", "0", 
        RowBox[{"-", "0.14277759330119488`"}], "0"},
       {"0", "0.9019078610942634`", "0", 
        RowBox[{"-", "0.2358817799647956`"}]},
       {
        RowBox[{"-", "0.14277759330119477`"}], "0", "0.07140307440990923`", 
        "0"},
       {"0", 
        RowBox[{"-", "0.23588177996485243`"}], "0", "0.06938393180670932`"}
      },
      GridBoxAlignment->{
       "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, 
        "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
      GridBoxSpacings->{"Columns" -> {
          Offset[0.27999999999999997`], {
           Offset[0.7]}, 
          Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
          Offset[0.2], {
           Offset[0.4]}, 
          Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
    Function[BoxForm`e$, 
     MatrixForm[BoxForm`e$]]]}], "}"}]], "Output",
 CellChangeTimes->{
  3.669853093636024*^9, 3.669853558352743*^9, 3.669853999483027*^9, 
   3.6698927975207243`*^9, 3.669895541676893*^9, 3.669896865789462*^9, 
   3.669899717678467*^9, {3.669900023163961*^9, 3.669900052250667*^9}, 
   3.669917636959354*^9, 3.6699185043723783`*^9, 3.669927444944366*^9, 
   3.669928978245573*^9, {3.669929019047792*^9, 3.669929043684908*^9}}]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Track a Falling Object", "Subchapter",
 CellChangeTimes->{{3.669896888318717*^9, 3.6698968927014217`*^9}}],

Cell[TextData[{
 "Repro an example from Zarchan and Musoff, ",
 StyleBox["Fundamentals of Kalman Filtering, A Practical Approach",
  FontSlant->"Italic"],
 ", Ch. 4, track a falling object, no air drag. There are three states: the \
position ",
 Cell[BoxData[
  FormBox["x", TraditionalForm]]],
 ", velocity ",
 Cell[BoxData[
  FormBox[
   OverscriptBox["x", "."], TraditionalForm]],
  FormatType->"TraditionalForm"],
 ", acceleration ",
 Cell[BoxData[
  FormBox[
   OverscriptBox["x", "\[DoubleDot]"], TraditionalForm]],
  FormatType->"TraditionalForm"],
 " ; and the observations ",
 Cell[BoxData[
  FormBox["z", TraditionalForm]]],
 " directly measure the position state."
}], "Text",
 CellChangeTimes->{{3.6698242225824003`*^9, 3.669824277939733*^9}, 
   3.669846748045307*^9, {3.669846787249598*^9, 3.669846816312376*^9}, {
   3.669846863222188*^9, 3.669846874990254*^9}, {3.669896945763245*^9, 
   3.669896947467445*^9}, {3.669917324348648*^9, 3.669917361292861*^9}, {
   3.6699174058075123`*^9, 3.669917417151168*^9}, {3.669917505829458*^9, 
   3.669917544340041*^9}, {3.6699274573924417`*^9, 3.669927506238433*^9}}],

Cell["State-space dynamics: ", "Text",
 CellChangeTimes->{{3.669824304050704*^9, 3.669824307618775*^9}, {
  3.669927512753031*^9, 3.669927515421135*^9}}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"ClearAll", "[", 
   RowBox[{
   "x", ",", "xd", ",", "xdd", ",", "\[CapitalPsi]", ",", "g", ",", "u", ",", 
    "F", ",", "\[CapitalPhi]", ",", "\[CapitalXi]c", ",", "\[CapitalXi]", ",",
     "G", ",", "\[CapitalGamma]", ",", "A"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"g", "=", 
   RowBox[{"-", "32.2"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"u", "[", "t_", "]"}], ":=", 
   RowBox[{"col", "[", 
    RowBox[{"{", "g", "}"}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"\[CapitalPsi]", "=", 
   RowBox[{"(", GridBox[{
      {"x"},
      {"xd"},
      {"xdd"}
     }], ")"}]}], ";", 
  RowBox[{"F", "=", 
   RowBox[{"(", GridBox[{
      {"0", "1", "0"},
      {"0", "0", "1"},
      {"0", "0", "0"}
     }], ")"}]}], ";", 
  RowBox[{"G", "=", 
   RowBox[{"(", GridBox[{
      {"0"},
      {"0"},
      {"1"}
     }], ")"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"\[CapitalPhi]", "[", "dt_", "]"}], ":=", 
   RowBox[{"(", GridBox[{
      {"1", "dt", 
       RowBox[{
        SuperscriptBox["dt", "2"], "/", "2"}]},
      {"0", "1", "dt"},
      {"0", "0", "1"}
     }], ")"}]}], ";", 
  RowBox[{
   RowBox[{"\[CapitalGamma]", "[", "dt_", "]"}], ":=", 
   RowBox[{"(", GridBox[{
      {
       RowBox[{
        SuperscriptBox["dt", "2"], "/", "2"}]},
      {"dt"},
      {"1"}
     }], ")"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"\[CapitalXi]c", "=", 
   RowBox[{"(", GridBox[{
      {"0", "0", "0"},
      {"0", "0", "0"},
      {"0", "0", "1"}
     }], ")"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"\[CapitalXi]", "[", "dt_", "]"}], ":=", 
   RowBox[{"(", GridBox[{
      {
       FractionBox[
        SuperscriptBox["dt", "5"], "20"], 
       FractionBox[
        SuperscriptBox["dt", "4"], "8"], 
       FractionBox[
        SuperscriptBox["dt", "3"], "6"]},
      {
       FractionBox[
        SuperscriptBox["dt", "4"], "8"], 
       FractionBox[
        SuperscriptBox["dt", "3"], "3"], 
       FractionBox[
        SuperscriptBox["dt", "2"], "2"]},
      {
       FractionBox[
        SuperscriptBox["dt", "3"], "3"], 
       FractionBox[
        SuperscriptBox["dt", "2"], "2"], "dt"}
     }], ")"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"A", "[", "t_", "]"}], ":=", 
   RowBox[{"(", GridBox[{
      {"1", "0", "0"}
     }], ")"}]}], ";"}]}], "Input",
 CellChangeTimes->{{3.669824321862136*^9, 3.669824412502771*^9}, {
   3.6698249550889482`*^9, 3.669824974432363*^9}, {3.669825095748258*^9, 
   3.669825198392029*^9}, {3.669825266512384*^9, 3.669825285148755*^9}, 
   3.669825329504478*^9, {3.669825374699018*^9, 3.669825375840673*^9}, {
   3.669825415034175*^9, 3.669825421974855*^9}, {3.669825460648127*^9, 
   3.669825478876505*^9}, {3.669825817207653*^9, 3.669825830270995*^9}, {
   3.669825930514666*^9, 3.669825965616769*^9}, {3.66982600325618*^9, 
   3.669826020647374*^9}, {3.669826054015787*^9, 3.669826181329473*^9}, {
   3.66982634781846*^9, 3.669826371993444*^9}, {3.669829154705411*^9, 
   3.6698292238043203`*^9}, {3.6698293241783*^9, 3.669829341615615*^9}, {
   3.669829382854457*^9, 3.669829406725457*^9}, {3.669829436839952*^9, 
   3.669829494105236*^9}, 3.6698296839787397`*^9, {3.669829811781673*^9, 
   3.669829875106605*^9}, {3.6698299323367767`*^9, 3.669829938959715*^9}, {
   3.669832105802023*^9, 3.669832186045053*^9}, {3.6698324031179028`*^9, 
   3.669832404197543*^9}, {3.669843236953289*^9, 3.6698432605816813`*^9}, {
   3.6698433266919193`*^9, 3.669843459166027*^9}, {3.669843640622575*^9, 
   3.669843675637025*^9}, {3.669849913082837*^9, 3.669849914111038*^9}, {
   3.669849980365006*^9, 3.6698500610568943`*^9}, {3.669850180383233*^9, 
   3.6698502056514473`*^9}, {3.669850236827403*^9, 3.66985030024732*^9}, {
   3.669850342469634*^9, 3.669850385651761*^9}, {3.669850599438808*^9, 
   3.6698507867904387`*^9}, {3.6698516557471113`*^9, 
   3.6698517501501703`*^9}, {3.669852452595969*^9, 3.6698524533392344`*^9}, {
   3.6698525124903727`*^9, 3.669852546383617*^9}, {3.669853059908752*^9, 
   3.669853067842987*^9}, {3.669854246353188*^9, 3.669854257450845*^9}, {
   3.6698542898916817`*^9, 3.669854293249909*^9}, {3.669891055201386*^9, 
   3.669891066142602*^9}, {3.669891240896646*^9, 3.669891254737685*^9}, {
   3.669897062487031*^9, 3.669897245431704*^9}, {3.669898846972335*^9, 
   3.669898888739764*^9}, {3.669899043676032*^9, 3.669899165229725*^9}, {
   3.669899249827806*^9, 3.669899325494771*^9}, {3.669899360503356*^9, 
   3.669899505986113*^9}, {3.669900106688882*^9, 3.669900127111576*^9}, {
   3.669901434589625*^9, 3.6699014885060043`*^9}, {3.669924922257389*^9, 
   3.669924935942034*^9}, {3.669925013626265*^9, 3.6699250277552977`*^9}, {
   3.669925069521391*^9, 3.669925076719947*^9}}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"ClearAll", "[", 
   RowBox[{"experiment", ",", "withTimes", ",", "myStyle"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"withTimes", "[", 
    RowBox[{"times_", ",", "obj_"}], "]"}], ":=", 
   RowBox[{"MapThread", "[", 
    RowBox[{"List", ",", 
     RowBox[{"{", 
      RowBox[{"times", ",", "obj"}], "}"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{
    RowBox[{"myStyle", "[", "s_String", "]"}], ":=", 
    RowBox[{"Style", "[", 
     RowBox[{"s", ",", 
      RowBox[{"FontFamily", "\[Rule]", "\"\<Candara\>\""}], ",", 
      RowBox[{"FontSize", "\[Rule]", "12"}]}], "]"}]}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"experiment", "[", "\[IndentingNewLine]", 
    RowBox[{
    "startTime_", ",", "\[IndentingNewLine]", "endTime_", ",", 
     "\[IndentingNewLine]", "timeIncrement_", ",", "\[IndentingNewLine]", 
     "aPrioriState_", ",", "\[IndentingNewLine]", "aPrioriCovariance_", ",", 
     "\[IndentingNewLine]", "observationNoiseCovariance_", ",", 
     "\[IndentingNewLine]", "fakeFn_", ",", "\[IndentingNewLine]", 
     "trueObservationFn_", ",", "\[IndentingNewLine]", 
     "trueObservationDerivativeFn_"}], "]"}], ":=", "\[IndentingNewLine]", 
   RowBox[{"Module", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
      "t", ",", "timesIter", ",", "times", ",", "fakes", ",", "ests"}], "}"}],
      ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"timesIter", "=", 
       RowBox[{"{", 
        RowBox[{"t", ",", "startTime", ",", "endTime", ",", "timeIncrement"}],
         "}"}]}], ";", "\[IndentingNewLine]", 
      RowBox[{"times", "=", 
       RowBox[{"Table", "[", 
        RowBox[{"t", ",", 
         RowBox[{"Evaluate", "@", "timesIter"}]}], "]"}]}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"fakes", "=", 
       RowBox[{"Table", "[", 
        RowBox[{
         RowBox[{"fakeFn", "[", 
          RowBox[{"timeIncrement", ",", "t"}], "]"}], ",", 
         RowBox[{"Evaluate", "@", "timesIter"}]}], "]"}]}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"ests", "=", 
       RowBox[{"FoldList", "[", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"kalman", "[", "observationNoiseCovariance", "]"}], ",", 
         "\[IndentingNewLine]", 
         RowBox[{"{", 
          RowBox[{"aPrioriState", ",", "aPrioriCovariance"}], "}"}], ",", 
         "\[IndentingNewLine]", "fakes"}], "]"}]}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"Grid", "@", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"{", 
          RowBox[{
           RowBox[{"ListLinePlot", "[", 
            RowBox[{
             RowBox[{"{", "\[IndentingNewLine]", 
              RowBox[{
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"Table", "[", 
                   RowBox[{
                    RowBox[{"trueObservationFn", "[", "t", "]"}], ",", 
                    RowBox[{"Evaluate", "@", "timesIter"}]}], "]"}], "-", 
                  RowBox[{"ests", "\[LeftDoubleBracket]", 
                   RowBox[{
                    RowBox[{"2", ";;"}], ",", "1", ",", "1", ",", "1"}], 
                   "\[RightDoubleBracket]"}]}]}], "]"}], ",", 
               "\[IndentingNewLine]", 
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"Sqrt", "@", 
                  RowBox[{"ests", "\[LeftDoubleBracket]", 
                   RowBox[{
                    RowBox[{"2", ";;"}], ",", "2", ",", "1", ",", "1"}], 
                   "\[RightDoubleBracket]"}]}]}], "]"}], ",", 
               "\[IndentingNewLine]", 
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"-", 
                  RowBox[{"Sqrt", "@", 
                   RowBox[{"ests", "\[LeftDoubleBracket]", 
                    RowBox[{
                    RowBox[{"2", ";;"}], ",", "2", ",", "1", ",", "1"}], 
                    "\[RightDoubleBracket]"}]}]}]}], "]"}]}], "}"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"Frame", "\[Rule]", "True"}], ",", "\[IndentingNewLine]", 
             RowBox[{"FrameLabel", "\[Rule]", 
              RowBox[{"{", 
               RowBox[{
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{
                  "myStyle", "@", "\"\<Position Residual / foot\>\""}], ",", 
                  "\"\<\>\""}], "}"}], ",", 
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"myStyle", "@", "\"\<time / sec\>\""}], ",", 
                  "\[IndentingNewLine]", 
                  RowBox[{
                  "myStyle", "@", "\"\<Position Residuals vs. Time\>\""}]}], 
                 "}"}]}], "}"}]}], ",", "\[IndentingNewLine]", 
             RowBox[{"GridLines", "\[Rule]", "Automatic"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"ImageSize", "\[Rule]", "Medium"}]}], "]"}], ",", 
           "\[IndentingNewLine]", 
           RowBox[{"ListLinePlot", "[", 
            RowBox[{
             RowBox[{"{", "\[IndentingNewLine]", 
              RowBox[{
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"Table", "[", 
                   RowBox[{
                    RowBox[{"trueObservationDerivativeFn", "[", "t", "]"}], 
                    ",", 
                    RowBox[{"Evaluate", "@", "timesIter"}]}], "]"}], "-", 
                  RowBox[{"ests", "\[LeftDoubleBracket]", 
                   RowBox[{
                    RowBox[{"2", ";;"}], ",", "1", ",", "2", ",", "1"}], 
                   "\[RightDoubleBracket]"}]}]}], "]"}], ",", 
               "\[IndentingNewLine]", 
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"Sqrt", "@", 
                  RowBox[{"ests", "\[LeftDoubleBracket]", 
                   RowBox[{
                    RowBox[{"2", ";;"}], ",", "2", ",", "2", ",", "2"}], 
                   "\[RightDoubleBracket]"}]}]}], "]"}], ",", 
               "\[IndentingNewLine]", 
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"-", 
                  RowBox[{"Sqrt", "@", 
                   RowBox[{"ests", "\[LeftDoubleBracket]", 
                    RowBox[{
                    RowBox[{"2", ";;"}], ",", "2", ",", "2", ",", "2"}], 
                    "\[RightDoubleBracket]"}]}]}]}], "]"}]}], "}"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"Frame", "\[Rule]", "True"}], ",", "\[IndentingNewLine]", 
             RowBox[{"FrameLabel", "\[Rule]", 
              RowBox[{"{", 
               RowBox[{
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{
                  "myStyle", "@", "\"\<Speed Residual / [foot/sec]\>\""}], 
                  ",", "\"\<\>\""}], "}"}], ",", 
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"myStyle", "@", "\"\<time / sec\>\""}], ",", 
                  "\[IndentingNewLine]", 
                  RowBox[{
                  "myStyle", "@", "\"\<Speed Residuals vs. Time\>\""}]}], 
                 "}"}]}], "}"}]}], ",", "\[IndentingNewLine]", 
             RowBox[{"GridLines", "\[Rule]", "Automatic"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"ImageSize", "\[Rule]", "Medium"}]}], "]"}]}], "}"}], 
         ",", "\[IndentingNewLine]", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{"ListLinePlot", "[", 
            RowBox[{
             RowBox[{"{", "\[IndentingNewLine]", 
              RowBox[{
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"Table", "[", 
                  RowBox[{
                   RowBox[{"trueObservationFn", "[", "t", "]"}], ",", 
                   RowBox[{"Evaluate", "@", "timesIter"}]}], "]"}]}], "]"}], 
               ",", "\[IndentingNewLine]", 
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"ests", "\[LeftDoubleBracket]", 
                  RowBox[{
                   RowBox[{"2", ";;"}], ",", "1", ",", "1", ",", "1"}], 
                  "\[RightDoubleBracket]"}]}], "]"}]}], "}"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"Frame", "\[Rule]", "True"}], ",", "\[IndentingNewLine]", 
             RowBox[{"FrameLabel", "\[Rule]", 
              RowBox[{"{", 
               RowBox[{
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"myStyle", "@", "\"\<Position / foot\>\""}], ",", 
                  "\"\<\>\""}], "}"}], ",", 
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"myStyle", "@", "\"\<time / sec\>\""}], ",", 
                  "\[IndentingNewLine]", 
                  RowBox[{
                  "myStyle", "@", 
                   "\"\<Position Truth and Estimates vs. Time\>\""}]}], 
                 "}"}]}], "}"}]}], ",", "\[IndentingNewLine]", 
             RowBox[{"GridLines", "\[Rule]", "Automatic"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"ImageSize", "\[Rule]", "Medium"}]}], "]"}], ",", 
           "\[IndentingNewLine]", 
           RowBox[{"ListLinePlot", "[", 
            RowBox[{
             RowBox[{"{", "\[IndentingNewLine]", 
              RowBox[{
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"Table", "[", 
                  RowBox[{
                   RowBox[{"trueObservationDerivativeFn", "[", "t", "]"}], 
                   ",", 
                   RowBox[{"Evaluate", "@", "timesIter"}]}], "]"}]}], "]"}], 
               ",", "\[IndentingNewLine]", 
               RowBox[{"withTimes", "[", 
                RowBox[{"times", ",", 
                 RowBox[{"ests", "\[LeftDoubleBracket]", 
                  RowBox[{
                   RowBox[{"2", ";;"}], ",", "1", ",", "2", ",", "1"}], 
                  "\[RightDoubleBracket]"}]}], "]"}]}], "}"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"Frame", "\[Rule]", "True"}], ",", "\[IndentingNewLine]", 
             RowBox[{"FrameLabel", "\[Rule]", 
              RowBox[{"{", 
               RowBox[{
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"myStyle", "@", "\"\<Speed / [foot/sec]\>\""}], 
                  ",", "\"\<\>\""}], "}"}], ",", 
                RowBox[{"{", "\[IndentingNewLine]", 
                 RowBox[{
                  RowBox[{"myStyle", "@", "\"\<time / sec\>\""}], ",", 
                  "\[IndentingNewLine]", 
                  RowBox[{
                  "myStyle", "@", 
                   "\"\<Speed Truth and Estimate\.7fs vs. Time\>\""}]}], 
                 "}"}]}], "}"}]}], ",", "\[IndentingNewLine]", 
             RowBox[{"GridLines", "\[Rule]", "Automatic"}], ",", 
             "\[IndentingNewLine]", 
             RowBox[{"ImageSize", "\[Rule]", "Medium"}]}], "]"}]}], "}"}]}], 
        "}"}]}]}]}], "]"}]}], ";"}]}], "Input",
 CellChangeTimes->{{3.669851430338047*^9, 3.669851502033433*^9}, {
   3.669899536592041*^9, 3.669899548535575*^9}, {3.6698996298270884`*^9, 
   3.669899633186829*^9}, {3.669900314852215*^9, 3.6699003162975388`*^9}, {
   3.6699007096190653`*^9, 3.669900729367638*^9}, {3.669900786936969*^9, 
   3.669900877284683*^9}, {3.6699011883365297`*^9, 3.6699012120691433`*^9}, {
   3.669901291581048*^9, 3.669901339352192*^9}, {3.669901538330483*^9, 
   3.669901561624837*^9}, {3.6699016052084427`*^9, 3.669901614758531*^9}, 
   3.6699160043532248`*^9, {3.669917585039975*^9, 3.6699176151858683`*^9}, {
   3.669917676738117*^9, 3.669917713022142*^9}, {3.669917758326228*^9, 
   3.669917804988469*^9}, {3.669917838306674*^9, 3.669917922370425*^9}, {
   3.669918394526479*^9, 3.669918496150769*^9}, {3.6699185263840237`*^9, 
   3.669918977085259*^9}, {3.669919052653118*^9, 3.669919244185959*^9}, {
   3.669920884342572*^9, 3.6699209622306337`*^9}, {3.66992145190014*^9, 
   3.669921575560782*^9}, {3.669921673850966*^9, 3.66992172296725*^9}, {
   3.669921759814693*^9, 3.669921831652063*^9}, {3.669922021269888*^9, 
   3.6699220481321383`*^9}, {3.6699222716218576`*^9, 3.669922383854618*^9}, {
   3.669922464237934*^9, 3.6699224796122513`*^9}, {3.669922640798594*^9, 
   3.6699227985424957`*^9}, {3.669922850086793*^9, 3.669922913662126*^9}, {
   3.669922944667591*^9, 3.669922987472431*^9}, {3.669923144226734*^9, 
   3.669923320184342*^9}, {3.66992385998766*^9, 3.669923914575399*^9}, {
   3.6699239506237297`*^9, 3.669923962503223*^9}, {3.6699239939831963`*^9, 
   3.6699240351519613`*^9}, {3.669924280308405*^9, 3.6699243849579897`*^9}, {
   3.669924662484598*^9, 3.669924668170545*^9}, {3.669924698491127*^9, 
   3.669924852660933*^9}, {3.669925129482019*^9, 3.669925257762621*^9}, 
   3.669925299497945*^9, {3.669925356043685*^9, 3.6699253631737337`*^9}, {
   3.669925422572912*^9, 3.669925434339156*^9}, {3.6699254867166433`*^9, 
   3.669925680370778*^9}, {3.6699257289241323`*^9, 3.669925767043036*^9}, 
   3.6699259990176287`*^9, {3.6699264139660187`*^9, 3.669926441139906*^9}, {
   3.669926497971624*^9, 3.669926588139385*^9}, {3.669926629925565*^9, 
   3.6699266547620068`*^9}, {3.669926727274036*^9, 3.6699267434169207`*^9}, {
   3.669926826613665*^9, 3.669926945976636*^9}, {3.669927071902417*^9, 
   3.669927105209777*^9}, {3.669929164920479*^9, 3.6699291912624683`*^9}}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"With", "[", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{"x0", "=", "400000"}], ",", 
     RowBox[{"v0", "=", 
      RowBox[{"-", "6000"}]}], ",", 
     RowBox[{"a0", "=", 
      RowBox[{"-", "32.2"}]}]}], "}"}], ",", "\[IndentingNewLine]", 
   RowBox[{"With", "[", 
    RowBox[{
     RowBox[{"{", "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"aPrioriState", "=", 
        RowBox[{"zero", "[", 
         RowBox[{"3", ",", "1"}], "]"}]}], 
       RowBox[{"(*", 
        RowBox[{"col", "[", 
         RowBox[{"{", 
          RowBox[{"x0", ",", "v0", ",", "a0"}], "}"}], "]"}], "*)"}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{"aPrioriCovariance", "=", 
        RowBox[{"1000000000000000", 
         RowBox[{"id", "[", "3", "]"}]}]}], ",", "\[IndentingNewLine]", 
       RowBox[{"\[CapitalZeta]", "=", 
        RowBox[{"col", "[", 
         RowBox[{"{", 
          SuperscriptBox["1000.0", "2"], "}"}], "]"}]}]}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"experiment", "[", "\[IndentingNewLine]", 
      RowBox[{
      "0", ",", "57.5", ",", "0.1", ",", "\[IndentingNewLine]", 
       "aPrioriState", ",", "aPrioriCovariance", ",", "\[IndentingNewLine]", 
       RowBox[{"1000000", "*", 
        RowBox[{"id", "[", "1", "]"}]}], ",", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{"dt", ",", "t"}], "}"}], "\[Function]", 
        "\[IndentingNewLine]", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"0", "*", 
           RowBox[{"(", GridBox[{
              {
               FractionBox[
                SuperscriptBox["dt", "5"], "20"], 
               FractionBox[
                SuperscriptBox["dt", "4"], "8"], 
               FractionBox[
                SuperscriptBox["dt", "3"], "6"]},
              {
               FractionBox[
                SuperscriptBox["dt", "4"], "8"], 
               FractionBox[
                SuperscriptBox["dt", "3"], "3"], 
               FractionBox[
                SuperscriptBox["dt", "2"], "2"]},
              {
               FractionBox[
                SuperscriptBox["dt", "3"], "3"], 
               FractionBox[
                SuperscriptBox["dt", "2"], "2"], "dt"}
             }], ")"}]}], ",", "\[IndentingNewLine]", 
          RowBox[{"(", GridBox[{
             {"1", "dt", 
              RowBox[{
               SuperscriptBox["dt", "2"], "/", "2"}]},
             {"0", "1", "dt"},
             {"0", "0", "1"}
            }], ")"}], ",", "\[IndentingNewLine]", 
          RowBox[{"(", GridBox[{
             {
              RowBox[{
               SuperscriptBox["dt", "2"], "/", "2"}]},
             {"dt"},
             {"1"}
            }], ")"}], ",", "\[IndentingNewLine]", 
          RowBox[{"zero", "[", "1", "]"}], ",", "\[IndentingNewLine]", 
          RowBox[{"(", GridBox[{
             {"1", "0", "0"}
            }], ")"}], ",", "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"col", "[", 
            RowBox[{"{", 
             RowBox[{"x0", "+", 
              RowBox[{"v0", " ", "t"}], "+", 
              FractionBox[
               RowBox[{"a0", " ", 
                SuperscriptBox["t", "2"]}], "2"]}], "}"}], "]"}], "+", 
           RowBox[{"gen", "[", "\[CapitalZeta]", "]"}]}]}], "}"}]}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{"t", "\[Function]", 
        RowBox[{"x0", "+", 
         RowBox[{"v0", " ", "t"}], "+", 
         RowBox[{"a0", " ", 
          RowBox[{
           SuperscriptBox["t", "2"], "/", "2"}]}]}]}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{"t", "\[Function]", 
        RowBox[{"v0", "+", 
         RowBox[{"a0", " ", "t"}]}]}]}], "\[IndentingNewLine]", "]"}]}], 
    "]"}]}], "]"}]], "Input",
 CellChangeTimes->{{3.669851430338047*^9, 3.669851502033433*^9}, {
   3.669899536592041*^9, 3.669899548535575*^9}, {3.6698996298270884`*^9, 
   3.669899633186829*^9}, {3.669900314852215*^9, 3.6699003162975388`*^9}, {
   3.6699007096190653`*^9, 3.669900729367638*^9}, {3.669900786936969*^9, 
   3.669900877284683*^9}, {3.6699011883365297`*^9, 3.6699012120691433`*^9}, {
   3.669901291581048*^9, 3.669901339352192*^9}, {3.669901538330483*^9, 
   3.669901561624837*^9}, {3.6699016052084427`*^9, 3.669901614758531*^9}, 
   3.6699160043532248`*^9, {3.669917585039975*^9, 3.6699176151858683`*^9}, {
   3.669917676738117*^9, 3.669917713022142*^9}, {3.669917758326228*^9, 
   3.669917804988469*^9}, {3.669917838306674*^9, 3.669917922370425*^9}, {
   3.669918394526479*^9, 3.669918496150769*^9}, {3.6699185263840237`*^9, 
   3.669918977085259*^9}, {3.669919052653118*^9, 3.669919244185959*^9}, {
   3.669920884342572*^9, 3.6699209622306337`*^9}, {3.66992145190014*^9, 
   3.669921575560782*^9}, {3.669921673850966*^9, 3.66992172296725*^9}, {
   3.669921759814693*^9, 3.669921831652063*^9}, {3.669922021269888*^9, 
   3.6699220481321383`*^9}, {3.6699222716218576`*^9, 3.669922383854618*^9}, {
   3.669922464237934*^9, 3.6699224796122513`*^9}, {3.669922640798594*^9, 
   3.6699227985424957`*^9}, {3.669922850086793*^9, 3.669922913662126*^9}, {
   3.669922944667591*^9, 3.669922987472431*^9}, {3.669923144226734*^9, 
   3.669923320184342*^9}, {3.66992385998766*^9, 3.669923914575399*^9}, {
   3.6699239506237297`*^9, 3.669923962503223*^9}, {3.6699239939831963`*^9, 
   3.6699240351519613`*^9}, {3.669924280308405*^9, 3.6699243849579897`*^9}, {
   3.669924662484598*^9, 3.669924668170545*^9}, {3.669924698491127*^9, 
   3.669924852660933*^9}, {3.669925129482019*^9, 3.669925257762621*^9}, 
   3.669925299497945*^9, {3.669925356043685*^9, 3.6699253631737337`*^9}, {
   3.669925422572912*^9, 3.669925434339156*^9}, {3.6699254867166433`*^9, 
   3.669925680370778*^9}, {3.6699257289241323`*^9, 3.669925767043036*^9}, 
   3.6699259990176287`*^9, {3.6699260750134974`*^9, 3.669926098888126*^9}, {
   3.669926148913107*^9, 3.669926163717383*^9}, {3.6699262040872517`*^9, 
   3.669926215419427*^9}, {3.669926972810597*^9, 3.669927047522114*^9}, {
   3.669927123898818*^9, 3.669927126154437*^9}, {3.669929215871078*^9, 
   3.669929246114154*^9}, {3.66992928796259*^9, 3.669929291625848*^9}, {
   3.669929341855577*^9, 3.6699293464151907`*^9}}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Thread", "::", "tdlen"}], "MessageName"], 
  RowBox[{
  ":", " "}], "\<\"Objects of unequal length in \
\[NoBreak]\\!\\(\\*RowBox[{RowBox[{\\\"{\\\", RowBox[{\\\"400000.`\\\", \\\",\
\\\", \\\"399399.839`\\\", \\\",\\\", \\\"398799.356`\\\", \\\",\\\", \
\\\"398198.551`\\\", \\\",\\\", \\\"397597.424`\\\", \\\",\\\", \
\\\"396995.975`\\\", \\\",\\\", \\\"396394.204`\\\", \\\",\\\", \
\\\"395792.111`\\\", \\\",\\\", \\\"395189.696`\\\", \\\",\\\", \
\\\"394586.959`\\\", \\\",\\\", \\\"393983.9`\\\", \\\",\\\", \\\"393380.519`\
\\\", \\\",\\\", \\\"392776.816`\\\", \\\",\\\", \\\"392172.791`\\\", \\\",\\\
\", \\\"391568.444`\\\", \\\",\\\", \\\"390963.775`\\\", \\\",\\\", \
RowBox[{\\\"\[LeftSkeleton]\\\", \\\"19\\\", \\\"\[RightSkeleton]\\\"}], \
\\\",\\\", \\\"378802.775`\\\", \\\",\\\", \\\"378191.344`\\\", \\\",\\\", \\\
\"377579.591`\\\", \\\",\\\", \\\"376967.516`\\\", \\\",\\\", \\\"376355.119`\
\\\", \\\",\\\", \\\"375742.4`\\\", \\\",\\\", \\\"375129.359`\\\", \
\\\",\\\", \\\"374515.996`\\\", \\\",\\\", \\\"373902.311`\\\", \\\",\\\", \\\
\"373288.304`\\\", \\\",\\\", \\\"372673.975`\\\", \\\",\\\", \\\"372059.324`\
\\\", \\\",\\\", \\\"371444.351`\\\", \\\",\\\", \\\"370829.056`\\\", \\\",\\\
\", \\\"370213.439`\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \
\\\"526\\\", \\\"\[RightSkeleton]\\\"}]}], \\\"}\\\"}], \\\"+\\\", RowBox[{\\\
\"{\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"1\\\", \
\\\"\[RightSkeleton]\\\"}], \\\"}\\\"}]}]\\)\[NoBreak] cannot be combined. \
\\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", ButtonStyle->\\\"Link\\\", \
ButtonFrame->None, ButtonData:>\\\"paclet:ref/message/Thread/tdlen\\\", \
ButtonNote -> \\\"Thread::tdlen\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.669929247950046*^9, 3.669929293281885*^9, 
  3.66992934752557*^9}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"MapThread", "::", "mptd"}], "MessageName"], 
  RowBox[{
  ":", " "}], "\<\"Object \[NoBreak]\\!\\(\\*RowBox[{RowBox[{\\\"{\\\", \
RowBox[{\\\"400000.`\\\", \\\",\\\", \\\"399399.839`\\\", \\\",\\\", \
\\\"398799.356`\\\", \\\",\\\", \\\"398198.551`\\\", \\\",\\\", \
\\\"397597.424`\\\", \\\",\\\", \\\"396995.975`\\\", \\\",\\\", \
\\\"396394.204`\\\", \\\",\\\", \\\"395792.111`\\\", \\\",\\\", \
\\\"395189.696`\\\", \\\",\\\", \\\"394586.959`\\\", \\\",\\\", \\\"393983.9`\
\\\", \\\",\\\", \\\"393380.519`\\\", \\\",\\\", \\\"392776.816`\\\", \\\",\\\
\", \\\"392172.791`\\\", \\\",\\\", \\\"391568.444`\\\", \\\",\\\", \
\\\"390963.775`\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"19\\\", \
\\\"\[RightSkeleton]\\\"}], \\\",\\\", \\\"378802.775`\\\", \\\",\\\", \
\\\"378191.344`\\\", \\\",\\\", \\\"377579.591`\\\", \\\",\\\", \
\\\"376967.516`\\\", \\\",\\\", \\\"376355.119`\\\", \\\",\\\", \\\"375742.4`\
\\\", \\\",\\\", \\\"375129.359`\\\", \\\",\\\", \\\"374515.996`\\\", \\\",\\\
\", \\\"373902.311`\\\", \\\",\\\", \\\"373288.304`\\\", \\\",\\\", \
\\\"372673.975`\\\", \\\",\\\", \\\"372059.324`\\\", \\\",\\\", \
\\\"371444.351`\\\", \\\",\\\", \\\"370829.056`\\\", \\\",\\\", \
\\\"370213.439`\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"526\\\", \
\\\"\[RightSkeleton]\\\"}]}], \\\"}\\\"}], \\\"+\\\", RowBox[{\\\"{\\\", \
RowBox[{\\\"\[LeftSkeleton]\\\", \\\"1\\\", \\\"\[RightSkeleton]\\\"}], \\\"}\
\\\"}]}]\\)\[NoBreak] at position {2, \
\[NoBreak]\\!\\(\\*RowBox[{\\\"2\\\"}]\\)\[NoBreak]} in \
\[NoBreak]\\!\\(\\*RowBox[{\\\"MapThread\\\", \\\"[\\\", \
RowBox[{\\\"List\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"1\\\", \
\\\"\[RightSkeleton]\\\"}]}], \\\"]\\\"}]\\)\[NoBreak] has only \[NoBreak]\\!\
\\(\\*RowBox[{\\\"0\\\"}]\\)\[NoBreak] of required \
\[NoBreak]\\!\\(\\*RowBox[{\\\"1\\\"}]\\)\[NoBreak] dimensions. \
\\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", ButtonStyle->\\\"Link\\\", \
ButtonFrame->None, ButtonData:>\\\"paclet:ref/message/MapThread/mptd\\\", \
ButtonNote -> \\\"MapThread::mptd\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.669929247950046*^9, 3.669929293281885*^9, 
  3.6699293475554256`*^9}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"Thread", "::", "tdlen"}], "MessageName"], 
  RowBox[{
  ":", " "}], "\<\"Objects of unequal length in \
\[NoBreak]\\!\\(\\*RowBox[{RowBox[{\\\"{\\\", RowBox[{\\\"400000.`\\\", \\\",\
\\\", \\\"399399.839`\\\", \\\",\\\", \\\"398799.356`\\\", \\\",\\\", \
\\\"398198.551`\\\", \\\",\\\", \\\"397597.424`\\\", \\\",\\\", \
\\\"396995.975`\\\", \\\",\\\", \\\"396394.204`\\\", \\\",\\\", \
\\\"395792.111`\\\", \\\",\\\", \\\"395189.696`\\\", \\\",\\\", \
\\\"394586.959`\\\", \\\",\\\", \\\"393983.9`\\\", \\\",\\\", \\\"393380.519`\
\\\", \\\",\\\", \\\"392776.816`\\\", \\\",\\\", \\\"392172.791`\\\", \\\",\\\
\", \\\"391568.444`\\\", \\\",\\\", \\\"390963.775`\\\", \\\",\\\", \
RowBox[{\\\"\[LeftSkeleton]\\\", \\\"19\\\", \\\"\[RightSkeleton]\\\"}], \
\\\",\\\", \\\"378802.775`\\\", \\\",\\\", \\\"378191.344`\\\", \\\",\\\", \\\
\"377579.591`\\\", \\\",\\\", \\\"376967.516`\\\", \\\",\\\", \\\"376355.119`\
\\\", \\\",\\\", \\\"375742.4`\\\", \\\",\\\", \\\"375129.359`\\\", \
\\\",\\\", \\\"374515.996`\\\", \\\",\\\", \\\"373902.311`\\\", \\\",\\\", \\\
\"373288.304`\\\", \\\",\\\", \\\"372673.975`\\\", \\\",\\\", \\\"372059.324`\
\\\", \\\",\\\", \\\"371444.351`\\\", \\\",\\\", \\\"370829.056`\\\", \\\",\\\
\", \\\"370213.439`\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \
\\\"526\\\", \\\"\[RightSkeleton]\\\"}]}], \\\"}\\\"}], \\\"+\\\", RowBox[{\\\
\"{\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"1\\\", \
\\\"\[RightSkeleton]\\\"}], \\\"}\\\"}]}]\\)\[NoBreak] cannot be combined. \
\\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", ButtonStyle->\\\"Link\\\", \
ButtonFrame->None, ButtonData:>\\\"paclet:ref/message/Thread/tdlen\\\", \
ButtonNote -> \\\"Thread::tdlen\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.669929247950046*^9, 3.669929293281885*^9, 
  3.6699293476131763`*^9}],

Cell[BoxData[
 RowBox[{
  StyleBox[
   RowBox[{"MapThread", "::", "mptd"}], "MessageName"], 
  RowBox[{
  ":", " "}], "\<\"Object \[NoBreak]\\!\\(\\*RowBox[{RowBox[{\\\"{\\\", \
RowBox[{\\\"400000.`\\\", \\\",\\\", \\\"399399.839`\\\", \\\",\\\", \
\\\"398799.356`\\\", \\\",\\\", \\\"398198.551`\\\", \\\",\\\", \
\\\"397597.424`\\\", \\\",\\\", \\\"396995.975`\\\", \\\",\\\", \
\\\"396394.204`\\\", \\\",\\\", \\\"395792.111`\\\", \\\",\\\", \
\\\"395189.696`\\\", \\\",\\\", \\\"394586.959`\\\", \\\",\\\", \\\"393983.9`\
\\\", \\\",\\\", \\\"393380.519`\\\", \\\",\\\", \\\"392776.816`\\\", \\\",\\\
\", \\\"392172.791`\\\", \\\",\\\", \\\"391568.444`\\\", \\\",\\\", \
\\\"390963.775`\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"19\\\", \
\\\"\[RightSkeleton]\\\"}], \\\",\\\", \\\"378802.775`\\\", \\\",\\\", \
\\\"378191.344`\\\", \\\",\\\", \\\"377579.591`\\\", \\\",\\\", \
\\\"376967.516`\\\", \\\",\\\", \\\"376355.119`\\\", \\\",\\\", \\\"375742.4`\
\\\", \\\",\\\", \\\"375129.359`\\\", \\\",\\\", \\\"374515.996`\\\", \\\",\\\
\", \\\"373902.311`\\\", \\\",\\\", \\\"373288.304`\\\", \\\",\\\", \
\\\"372673.975`\\\", \\\",\\\", \\\"372059.324`\\\", \\\",\\\", \
\\\"371444.351`\\\", \\\",\\\", \\\"370829.056`\\\", \\\",\\\", \
\\\"370213.439`\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"526\\\", \
\\\"\[RightSkeleton]\\\"}]}], \\\"}\\\"}], \\\"+\\\", RowBox[{\\\"{\\\", \
RowBox[{\\\"\[LeftSkeleton]\\\", \\\"1\\\", \\\"\[RightSkeleton]\\\"}], \\\"}\
\\\"}]}]\\)\[NoBreak] at position {2, \
\[NoBreak]\\!\\(\\*RowBox[{\\\"2\\\"}]\\)\[NoBreak]} in \
\[NoBreak]\\!\\(\\*RowBox[{\\\"MapThread\\\", \\\"[\\\", \
RowBox[{\\\"List\\\", \\\",\\\", RowBox[{\\\"\[LeftSkeleton]\\\", \\\"1\\\", \
\\\"\[RightSkeleton]\\\"}]}], \\\"]\\\"}]\\)\[NoBreak] has only \[NoBreak]\\!\
\\(\\*RowBox[{\\\"0\\\"}]\\)\[NoBreak] of required \
\[NoBreak]\\!\\(\\*RowBox[{\\\"1\\\"}]\\)\[NoBreak] dimensions. \
\\!\\(\\*ButtonBox[\\\"\[RightSkeleton]\\\", ButtonStyle->\\\"Link\\\", \
ButtonFrame->None, ButtonData:>\\\"paclet:ref/message/MapThread/mptd\\\", \
ButtonNote -> \\\"MapThread::mptd\\\"]\\)\"\>"}]], "Message", "MSG",
 CellChangeTimes->{3.669929247950046*^9, 3.669929293281885*^9, 
  3.6699293476599693`*^9}],

Cell[BoxData[
 TagBox[GridBox[{
    {
     RowBox[{
      GraphicsBox[{{}, {}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0., 0},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Position Residual / foot\"", FontFamily -> "Candara", FontSize -> 
            12, StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Position Residuals vs. Time\"", FontFamily -> "Candara", 
            FontSize -> 12, StripOnInput -> False], TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0., 1}, {0, 1}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.02], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}], " ", 
      GraphicsBox[{{}, {{}, {}, 
         {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
          0.008333333333333333], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1: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
I/PF8/pxzasyXzyvLwdNkPnief25f5TMl6DOH3aZL4HXr9d1Ml8Cr2/fSJD5
Enj92z5Y5kvg9XF+qMzHvKCeATJfAq+vN3vIfAm8/h7cSeZLIKchywJeue+N
fAnU54UpabMveSNfIinXLTzijXyJ5Pveod+X/OCNfIl0f6Eu4YNCb+RLpPG/
2Dw/XAJebiJNSbv2zIezwMtL5PbIfwG8wkSKdfuu91ITeMWJ1PJ85TeLYsAr
S6RRld98NreXzJdI5t31bVNdZL4kzpdT54V8SXQnItbNctkL+ZLos6Ohu+KO
eSFfEgnc0V5FXsiXRId+32xy2uiFfEmk3FW9DLzCJFLiH5sNXnESXTrV/PyW
l8ArSyJPIfBdK3jVSfTngjsRL2rB0yQzTxsi8yXTQPHCru4yXzLtnD1tVXl9
N+RLptPNz1duu9oN+ZJJ6N4870Q35Evm99WVdEO+ZPqu9wtTXDd1Q75kCt01
e1rZCvDKkkk0b/nKueBVJ9Mfm027R0wAT5PC7eFpBy8ohX9fGg+eNoUGi/sX
9ZX5Uqif4Md3lflS+Pm3GzyRL4XzfXrNE/lSOL/pd0/kSyGhr/JhiSfypZDA
rSrc5Il8KexH9wF4mlT2VzkXvKBU9rtkAnjaVPIQ/sPs4DlSqVw84EA8eLmp
3H5j+8p8Ku+Bu8yXyu3/boMH8qXSTdE//K55IF8qbRD958sTHsiXyv0rusQD
+dL4ffd87oF8adw/dSvA06bR80qHmwOeI437t3k8eLlp3P9LreDlpXF7GLTg
Fabx+NnfR+ZLo5fE+Epwk/nS6IEYfzsedUW+NB6fva90Rb50zldwrCvypfP4
fvxjV+RLp07K+N/QFfnS6T+6hPEnl4KXm04ffjrTV/smeHnp7G/tv8ArTKdZ
4oFNJvCK0ykzYfwvo2LAK0unLcoADJb50ulL0+56ty4yXwbzJtS6I18GiXjX
isrdkS+D/EWDdj7qjnwZNG3V/YXP/eCOfBk0UgT8fJ078mXw+95+F7zCDBLN
/+mw18ErzqClAeWnZo8FryyDRLrZJXrwqjNIvO579UPA02Rye8QEyHyZ1H5h
wZ1XNTJfJn0kCuDmGjfkyyTRPE4XL7ghXya9pXSYg27Il8n5kra7IV8mie5y
YcpqN+TLpHhRQP/3P+CVZZJdFJifp4JXnUkNbRcWXH4WPI2O/XWkgRekI9Gd
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

           "]]}, 
         {RGBColor[0.880722, 0.611041, 0.142051], PointSize[
          0.008333333333333333], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1:eJw11nlYjPv/x/GU9tKmPVq0l7Zp3+Y1NU3TzHSKQh2yHXv2rMdOJ3JCOgdF
ItuxRAjJ4ShEKEJRkTUJrdq1/fr+en/uP3Ld19z3437dz2uuuZjOXDJ+trSU
lBQG//zvX3Zs9jjx6d5iPd7BtP8deVzeiOevrIPZ+UNuzxWVGb946vFcOf87
XnBLnfbbhzuxzyu54i/i+c8Ferwh7T13nZ2O9SJvdv0nrtqG523XA/R4c+f8
76jl5sSenR1vz+7/ynXiHHsZYKnHKyn+31HP3Tc3OvulAfOaue82/LbxwjB2
fQt3IFlWurlHl/xW7uyQA0uihrP727ge6/K1D3Tp0vM6uKsbF3mbd+qS18Xd
r7kn5XizLnnd3FGjPdrnVDCvhxt0OTJG66Eueb3c/Ze2Rxj/x7x+7qHb28dP
KGCeFE6rZfo9uc2ul4Ky8bCL+ueYPwz6BvV//3GK3S+NM1IuOgbp7Hky6Ex6
LYzcx7zh0Eh9png5kXnDkb3L6Gn9NubJInJ8mXbKIubJIXv31C2xU5gnj2PS
82oVxzFPAWdNGt3GcJmngOnDJfOyOMxThOIjhZD2scxTgsHpI/s2WDBPGd96
9o5fbsQ8Fey5XKybMpJ5Krig9XL1a1XmqUJN48k/YxWZNwKNH59dDpNlnhpC
wqR0PnfqkKeO4x46hte/69D16li93jEp+b0O+epotlxxsK5ch/f/t3M0sKg2
8+vxBzr0PA082mK6oquAeZrIOJbnZXWNeZp4ahV6Mv4a8zTx82ZJ5rPzOrRH
C/vrLk3hnGSeFoKKta81pDJvJFoMK/ztU5g3EgnHi7ZZxTNvJDZfizSK/p15
2tB64rm2M4552vCTfLPKms08Hdy/P+1AVgzzdFDZv6h5/gTm6cBhVIqBp5h5
ulg/acEGBR7zdPHVNtrFk8M8PTSurNpw0YZ5egibrPamxZh5euA8XZy7R5d5
+pjCWdappcw8fTx+XnvWt1+bPANcsdVsTm/VJs8AX7oyuVe/apNnANf2sRr+
NdrkGcJq5qpJS6u0yTOEzLD7JXtLmWcE5asqNU5FzDPCnZWiWbU3mWeEW1Gn
vztcZd4omMa1iIPOMm8U9p5W+Tf1CPNG48D8H6/1DjBvNGY6Bj802M280TAS
dCp7xjPPGAklnOas9cwzhnBn3b+pS5lngmT9uaVP5mjT98sE8Vey4378ynwT
7JOpaXaJYPeboGuxlUVoCHueCeY2zPq30I95piitkDKsdGWeKbwm3lo+xZZ5
pnA86Fu53Ih5puh4n62Vq8E8U4xfsTH5ohzzzHDfQk/qfO9I8sxwccz1ncY/
RpJnhtac60l760aSZ4bJSZ/93lePJM8Mym/kXNaUjSRvDMRWj/jdj0dSnzEI
j5zyfnwB88YgwnilzZJc5o3BgxFHHdZcYN4YWNsWqR89wTxzbPkv7lXhIeaZ
o+Z+ab75X8wzBz/p0WvbROaZQxy6bsLyjcwzR3PNHzoT45hngSKZCDgvYJ4F
3pVJbtTMYJ4FFO7ERrVGMc8CQjczlyNhzLNAn1+q2Wg+8yxxbmXo7w3ezLNE
WuFfZnLOzLPEtdCzz/WsmWeJuP0HUotHMc8SD5fInXisyTwruMu7JNxUYp4V
bmsPX3RDmnlWqJpS8yKxW4s8K3zeEb6vqkmLPCtkVWdW7qzVIs8aoaVjcw5W
a5FnjeFh2Y6lZVrkWcNtTOOIN4+ZZ40NWY1TjxQwzxoBCf7R3bnMs8EU28Df
zl9gng1K8u8oJJ9kng1cIkWLog4xzwaHQ3d9fpHCPBvk7nL7vCCRebawdJ2h
Xr2Rebb4fZPBtYMrmGeLNplbJYkLmGcLa/+q8/unMc8WGRLd9zWRzLNDwJvd
vAUS5tmhy/7J2I9gnh3GWXyI1/Jknh2qy5f/9vdY5tkhYm7ltDoz5tkjf5/y
wnOGzLPHlZr135dqMM8eDoeb4l/IM88eCxtWC2b1apJnj7ucDJl5LZrkjUV7
dWzF7DpN8sZCPfPO9B9vNckbC80Yu72cck3yxmLclIelno+ZNxYHc4x+Lipg
ngNUlp7TGbjGPAf88O1U6MtingP6jxWen57JPAfoZESa3jvAPAfI2caIgnYx
zxFeM2Sa07YxzxEOez52Kq9lniPWzArv11rCPEfkLJP1XD+TeY4Ijdt8IWsS
85xQCMPnA6F0vb4TLPvPnIjha9LvjROksy4Om+5N90uccPPP/cIjjux5TlBf
MmHrB3PyNjth3idXm1QD9nwnuF3W6xlQIy/HCd99V2uPk2N7nCAeJR2Q06sx
5NU6wfek96ZTPzRonzOcYqf9F1inQfucIbc2qHfqWw3a54ysMsNpBWUatM8Z
ZfNMt0x/RN4cZ5SMc1zmn0/eZmcY7rD8RymXvDRnZPQbPRZlkZfjjFq3b22f
Mskrdsb7PXFqOw+wfc44r8eXzEli+1wwoaQhJm4L2+eC0sbyvuxVbJ8LEo8U
FlxYyPa5ILDq3iG9mWyfCw7d6U4xnsT2uaAu9x/rWAnb54LO8pRk40C2zwWf
Fn+xme3F9rkgwOReqq8j2+cCbsKcspmWbB8HblGh344Ysn0c9DyJuDhBQ4O+
PxzkHLrekSDP9nFgnq5R2dirTvs4aO2wM3v7Q532cWArSbeLr1OnfRwcijg/
aWO1+pCXM+iJX97697k67eOg7v6ylXcekFc7eH5jlvjqLfKkXOHz12nrX3PI
03eFgeWlyTNPk8dxhVaIhuvNDPIkrvj49JHTv3+zfa7Im2X6YeVOts8VY3MP
z+3YyPa5YpezlVruCrbPFVs+Xtrybj7b54rzMQ+LxTPYPld0+6x0nDWR7XOD
6u5mX00J2+eG8XIG/VPB9rnBWnVWbaQ72+cGeR4qNOzZPjfYKKr0CM3YPjdE
faoW8HXYPjfY2/9MDlNm+wY/v2SaPX4Y2+eGRWdPeEU0qNE+N8Q+W6WfXKlG
+9zR3RzeeeqeGu1zhwQWF1rPq9E+d9wNefW4IlWN9rlD5zPPQTmBvDnueDBr
QfiNJeRtdsfTsTLHp/9KXpo7bsgv/ZwXRF6OO/6LT46VOJFX7I7XpT/ybxux
fe54ZzJn7Gl5ts8DAzfjbwl+jKB9HvivYNr7zDcjaJ8HtGcpTql4MIL2eSD3
OO/xsysjaJ8HTPJrfxUeGUH7PJDisGP58T/JS/OA/UZJ+6IV5OV4YLlu097m
GPKKPfCHoXWEi5C8Wg+MypiANy7kSXlCM8nJNsOY7fNExUTLNXeV2D5PzGqa
GLyyU5X2eSJetmn7vE+qtM8TfUVb+q+XqNI+TwSMGSdYf02V9nmi/JRJgU2m
Ku3zhPqBxhU1O8kr9kSTbUfCrhXk1Xoiqi/q+9up5El5oSO+19pQRJ6+F850
Chf3ccjjeGFqisv0jcZsnxdWnbHjLFBi+7yw/5PVhjHtKrTPC+lxfh+uvFWh
fV5IELzMn1msQvsGP9+n/nfhVRXa54U9Dx3ruzJUaJ8XCud7OwQlkifljS7f
F9u0l5On743MA833BiaTx/HGnNc8zRFC8iTeEPzn6dPGIW+ONzouTXqlbsz2
eaOposCmWJHt88bBc/c7tNqVaZ83Fhhcr3n7Tpn2eeP5Z9PRD4qUaZ83pl2O
/hByWZn2+eBqj0la+SFl2ucDfYHH87IE8jg+0Lv6LsghjjyJD2b8ujpMdip5
c3xwManmVU0QeZt9sHSE7NxjjuSl+WDZKf/0+YZsnw+mxq7tOyzL9vngdObW
E9NblGifD95OdfRvfq1E+3zx9efVsn2FSrTPF3e7ApfdzVaifb4YbnZn+fc0
Jdrni+wlFnKZf5A3xxfrfh3DG72EvM2++GPnUdfEyeSl+UK06eaW8CDycnzR
eVFhe6cjecW+6D15Y6mJIds3+PmZFO1kWbbPD/duzT+2rkWR9vlh9lbtX2+8
UaR9frh1KLVp6X1F2ueHBxXzs8wuKtI+P8jtiq7ITFOkfX54l+lbbB9PXpof
ooRTM28tJi/HD4tvhf/aHU1esR8MdqQoRfHJq/XDhdBwHRVH8qT8EfFRtbpc
j+3zxzrXctsFMmyfP5wnL+z6VK9A+/zRFuvnoPpKgfb5Q2v0wUjZAgXa549r
pXVik7MKtM8fTU6Pd//ylwLt88dvumYKm9aTV+yP7MOqj6rmklfrD43ySX5H
wsmT4mJDpNuEuT7k6XNxv2Rr7GUL8jhcSKr0EvapsX1cjIsILqrqkqd9XARO
H1h/+aM87eOi7vimxjXF8rSPO/j9EiysvSZP+7hwzVh3x/yoPO3jYo2U3Izi
RPJquWhMPv/DbyV5UkDnm8brC6bR9SqAabaJZpiQfH0gQ/BLdLAL3W8JRHtI
LGSN5On/N8BCVelGbznyAEzRUsg52iJH7wOsvltRx38tN+RFA3+Oyzc9UChH
7wfUhu74O/+C3JAXB4zya163LE2O3hfYcKax0C2evCTg+rSqGQsXkZcG5FU/
3ts1ibxTwJzbPy/LBpCXA1TckF2WaEdePrBr3bkVXdrkFQMrT9jLHRuQHfIq
Ad4x3897vspSL2C2xS8a8s9lh7xWYO6mKae5N2WpHw8iLbmUplOy1I+H6PCD
Jycmk6fPw3dB3s5Lv5NnyYOG3jBHq9nkcXjoPnTg8r5w8sDDP01uQQt9yJPw
EN6392SMOXnRPOytHSloUCVvDg/pulsVbnQMp348zIsWBWi+H079eCjZVHPk
9cPh1I+Hxcf+cUzIGU79eCjfV195JX049eNhSnf//osJ5OXwYH1B2nvdEvLy
eeC1ZQdERJNXzENlkdzn7gDyKgffP9241c2evFoehKudUyu0yWvloWhVTpxo
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
WOR9kHeC+knQYpf4h3p+Ov4PNhuVbg==
           "]]}}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0, 0},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Position / foot\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Position Truth and Estimate vs. Time\"", FontFamily -> 
            "Candara", FontSize -> 12, StripOnInput -> False], 
           TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0, 57.5}, {0, 400305.22230612487`}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.02], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}]}], 
     RowBox[{
      GraphicsBox[{{}, {{}, {}, 
         {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
          0.008333333333333333], AbsoluteThickness[1.6], 
          LineBox[{{2.0925117552437476`, -64.93470891795505}, {2.1, 
           53.2011604253139}, {2.1239183758862525`, -64.93470891795505}}], 
          LineBox[{{11.720017351984549`, 66.67060376374002}, {11.8, 
           60.600055741977485`}, {11.9, 45.11818023131127}, {12., 
           22.023651913079448`}, {
           12.100000000000001`, -3.6408049394904083`}, {
           12.200000000000001`, -25.312925509376328`}, {
           12.3, -19.897277385566667`}, {12.4, -26.82946683544924}, {
           12.5, -40.85332571355866}, {
           12.600000000000001`, -58.289692494468}, {
           12.700000000000001`, -53.228587424801844`}, {
           12.740332459273313`, -64.93470891795505}}], 
          LineBox[{{13.246225622302738`, -64.93470891795505}, {
           13.3, -54.93230737243357}, {13.4, -41.7416007055308}, {
           13.479743321216235`, -64.93470891795505}}], 
          LineBox[{{13.520466736972445`, -64.93470891795505}, {
           13.600000000000001`, -42.04010847735026}, {
           13.700000000000001`, -64.76451675210501}, {
           13.8, -55.80492447890265}, {13.9, -43.93098898050903}, {
           13.941132434071397`, -64.93470891795505}}], 
          LineBox[{{14.344229711593782`, -64.93470891795505}, {
           14.4, -55.53056161463792}, {
           14.451128952393413`, -64.93470891795505}}], 
          LineBox[{{15.191024915463137`, -64.93470891795505}, {
           15.200000000000001`, -62.5936914265767}, {
           15.229605545650674`, -64.93470891795505}}], 
          LineBox[{{15.381495379627294`, -64.93470891795505}, {
           15.4, -63.670795779441505`}, {
           15.418064075646752`, -64.93470891795505}}], 
          LineBox[{{15.575894720740777`, -64.93470891795505}, {
           15.600000000000001`, -63.113849392050724`}, {
           15.682435692452229`, -64.93470891795505}}], LineBox[CompressedData["
1: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
suDH64WlZ/b8NNqE6ufPMjMbz0d7q1D94Bm9JoFrhFQ/2GF/4cNshuoHS37o
fdLSWMHXj7t++Jm0yoyffyw8dCHT+IBiPcEr6vPLf7Lg52MUd/1+k6p5ZjR/
4Z4bs8ot1PjzYQosivzfill91O/h0Nh7V0f18/u1FD49orc+uJf2dy7PxvqK
D938/RbB7JbGL/E1fL5zsCw98aNYg+oHF8cUt29SpvkFx0ZqaFuoUP1gleMh
x9yov3dx968uC14ppPrBqUP9Y5wNqH4B+L3c7kr7S/55CGCJf4eytI32O1j2
T9Tmn5/x/UkMZ2bclasb0PyDHdkP0W3U/1g4eb1hiUSf5h9cOTOs6awVzT/u
+2xNzjw8TPsjN75v5G9XnvP3mwLfNjfel0j7rwTONgoz8e3n80m5fN8JvrNV
9Ffu+9pab95T4Z9HESwqDlldqUzndW78+Zf/WUfPQwbXBHTVHlWi/vbv980e
Okb9ppnLJxFvy5XS/OO+LzuvwbuTz9fL/d2j4eoZPep/0zHfzQbHNdF+KIBf
Vt8vvaNG/RAWReffH91G9YOZRKtHhwbofAgnl1a6xBjR/OOu1z/w9NUYmn9w
Dduck6ZP8w92iyyNNDOh+QezLzZ6hOnx/TYFlhju63jWTPOPy3ezRlUs4MeT
wqba/6mVaVF/5b7PzDIvVpPqx+U72ff3Ix3aD2Dhgnnh+Wp8Phl3fezSoT80
qX5wblmEl1iZ6gdnR1w/OPUj7Wewd9SbjBE6P/TCBVlH33oN8/mYGchXmKb6
cZBfHwKYedYZ+LSS5h/nxRf7K3upfvDMyW/KdGg8D5ht6Qg/84Ufj+XGa/rz
RuwwP14ILIu7dvbIX/x4UbCbY71boRrVD859UG7vqU/1g4vdr+uEa9F+BNtW
DgdWGVH9ZnDzWzfPV4fqB8uVi37LHU31g3uC77z9a6iCrx+X/89j2q23+POV
jMtrLI1c+Y7O77Dk9ULfn3r5+dwMh540+9F9LK1feFl2mF6GmN9feuGWPZpO
b51p/QayTNz7R0U7nKj/wQ4dE1P97an/wSPfBab10/4jhuvSeztm03nJA64U
xTuEKdYvvPf8koXLx9P6hQ+FH/bXcKX1CyebZGSHudP+wf39zpdpU1z58VLg
FU0xgdsdaf3Cu4qiTRpoPCmc5q/c0m1H+wfnyq0bauxp/4CLBiKlbXS/5+Cs
rNWdmXT+k8Gl+fN78hX7K6y2+0JmrAP1Py6f6/mhM+OofnCkand1iYjqB384
nmPgbEf1C2KZLUu3rqyj87EATs2VxZlTvzKDpwatbBi0ovrBQ57C8PVjqH7w
C1ct5osp1Q8+HeBfN8BQ/eDQdV1Juv38842CGV2Dt4FldL6F2eTHSo1dfD9I
gUUCj8b73bR+uetrRplbnKT9A37ZuOq5rybVDxY4xMTX0Pm7CFYpsL0yzpjq
B0cI4qfn034pg/v0rFiJCdUPNt1UuPugAdUP3lFpsCxLcZ6Hs+PKrt6k++2F
NTJ9/XYq3j+C8ds1IXY07ecC2Ltx8cvjDNUPzs4//aaWzvNiWFJu2m5M5zkP
uOfvPdcch/n1wcKiqY7jHTvofYq7fvPzOUm0PqLg2IJH8nQ51Q+OVDvqH03n
uRR4RTZz6rU+zT84rcRrfd5omn+wXHVF6xw6D+bBbNuevWe1qX5c3nsqW7SV
qX5cvgMa4ppPtH5htzftG+qH6fzC3X/A/ulBd/n+0szlXTet9Ndq2j+48eMn
Dp94z99vLywL+eazTR29/87E57urNjQ+pfrBofuc1IWf+fs1g70v3HDLoPOV
mHOv3c6rajT/YFGTy4buHqofLCwvfVLzkeoHy46Y/vF9A9VvJrcfvCu8/I7q
B0sGEw9vvE3zj/t8kGaRoTr1P7i4Vf9zpzL1P7hgW8Kb+pd0fuG+X/Z939sh
vj8XwRrymPg+Ldo/4LL6nCp/6qcyOE3bbb2P4v0RDpk8lH5FSPsHLDdzFF9W
/L8FZr+9+YvqU76f9nLfVyg0i/+H9o9Z8JzB/S0faP+AWaMX8RPaaf+Yxb1v
uYUFf6D5x/19/5vdBr00/2BJ7WjP5naqH9zD/Ch3kVP9YCZn62dDmi9R8OJp
25JU9an/cZ9f0uKaZEjzD069NXvERTH/YLcThg/DjGj+cX8/WXqmXHF+hiP8
Jcr1NF4RbDvPee44PZp/3PXiECadzmsyONP/jGwbjVcN92Wtv2xtTOuXyzui
I7an8bpgeYP9kJKA1i/cdirY9IE6rd/Z6Jf+841UaX0I4OqdA1raRrR+4Rq3
mMoseh8Sw5dV1xXdsKb5B6t4qjo/p/cjFhYOJ53Ps6T+x413PSztFv0/NQqW
i/wdm8ZQ/eDK226abpZUP7hOLg9Op/c3Cfzibv4HP3OqH5dn6bXru62pftz3
v/5m4kkarwhevcpGrdyE6gd7mmSE/2NJ9YO7HOJ7NcZR/eCjT6Ze0KF8zbD3
tHdp3wqoftz9PGaeNehQ/bj8W19nF6tQ/eag/2yK+fkrIdUPdoho+tVFmeoH
p1Wlt0xTp/rBKrmLjlfpU/1g0xJ7A08B1Q9O2NVa8lSxf8BH64K2tOtR/eC6
mmbLbYr3NzjTp79UZET1g08XNqYMK+YfbDFrsbe9AdUPZupsc3QV/Q+OfdS7
9Yyi/8GhmtvHHlL0PzjZ6Gnywj7qf7Bs055rCzvo/ALLU3Y+OaBC6xfO9tpg
Fa5K6xcuKzr7+xglOv/BGnlbj7mo8v2ACcF687rTcpv6iwBuu/REQ2OI1i8c
e9ovLrubzn+w5PON/76+TuuX80a92XE36PzHfT7z/OkEOk+GwMKEht2tSnR+
5q5/tUk0f4TeP2CNQfMrX1To/AfLX5Q4lWhQ/+M+7/hf7RF16n9wzYUGD2u6
3zx4pmnqhlG6dP7j8s45NfalOvU/uOBbk9TRA/T+AUfatYS915H9P7pWJuE=

           "]], 
          LineBox[{{0.06980724539530017, -64.93470891795505}, {
           0.07000798214095084, 66.67060376374002}}], 
          LineBox[{{0.18141648239185773`, 66.67060376374002}, {
           0.1819614033703152, -64.93470891795505}}], 
          LineBox[{{0.9714747224919169, -64.93470891795505}, {
           0.9742389118603176, 66.67060376374002}}], 
          LineBox[{{1.1072688660423153`, 66.67060376374002}, {
           1.1109200277687687`, -64.93470891795505}}], 
          LineBox[{{1.6100106840556176`, -64.93470891795505}, {
           1.6233854824002296`, 66.67060376374002}}], 
          LineBox[{{1.8556191765092294`, 66.67060376374002}, {
           1.8633433875281866`, -64.93470891795505}}], 
          LineBox[{{2.6271844482673155`, -64.93470891795505}, {
           2.637445384416479, 66.67060376374002}}]}}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0, 0},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Speed Residual / [foot/sec]\"", FontFamily -> "Candara", 
            FontSize -> 12, StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Speed Residuals vs. Time\"", FontFamily -> "Candara", FontSize -> 
            12, StripOnInput -> False], TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0, 57.5}, {-64.93470891795505, 66.67060376374002}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.05], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}], " ", 
      GraphicsBox[{{}, {{}, {}, 
         {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
          0.006944444444444445], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1:eJxV13uMzXcexvGzlmXtWKqjOuo2Y+5mzO3Mfeb8fud+jYaVihCNrSWNxiVj
bZXICBVqokQrZoUQa0KE6E7KmDGTCal2lcZWukorVVK6VGo7Vday3fn93t/n
j50/iHHOc57zfLwcMn+/ZPofBng8Hrv/B+fn//taeqL3z63O10lrwpt3p9U0
6dd/sxqHXt6/dMWJXm+F83XJmvP26QcHV+r3r1ir0o9Gv159opeg65bz/Yxm
Pf6m1TG+P3Ddid6FC5yvW1Z/WuNbG/T8f1o/5fcnbjrRe+G88/Wd5cQ9blHe
fcuN2arH/8ua3h+4aLvy+6xlTsEdev6PVn+9/C9b9Xo/Wf1xb6TvVt4jy3lc
aq/y/m19N60/cL/y/mM5b7e7TXlPrALnDR9S3n+tWH/clCPK89jut4/p8R67
/9127GlX/i/svyztL3hczx9gn3mQP2d4p17vl3Z/3NFot/IG2s6Pzb3KG2g7
9+g4o7xBtvN2759V3q9sp17+OeUNtp17zLugvCG2+7yLyhtiO/f4+yXl/dp2
7jH0svKG2s49AleV9xvbuccb15SXZjvfbb+uvDTbucfdm8obZjv3mHRbeb+1
nXpz7ihvuO3EvXNPeSNsN+e+Hj/Cdu4x6Eflj7CdezQ+PNHrPr3iGdu5xx8f
6/WesZ17HH2qvJHuPW55OkzeSPce4wd2mLyR7j1mDu4wfZ517/H20A6T96x7
jw+HdZi8dPcenmeUl+7eoyZdeenuPZaOVt4o9x4HxyhvlHuPr8cp7zn3HhmZ
ynvOvce0bOU9597jrTzljXbvcbpQeaPdezwuVt7z7j0qypT3vHuPRV7lPe/e
Y3+18jLce3xZp7wM9x7pPuWNce+R8itvjHuPN0PKG+PeozuqvBfcezxIKO8F
9x5TpipvLD6mKW8sPmYobyw+ZipvHD5mK28cPuYqbzw+5ilvPD7mK288PhYq
bwI+FilvAj4WK28iPpZ1mD9fE/GxXPkT8fG6nj8RH6v0ehPxsUZ5mfhYq7xM
fKxXXiY+NiovEx+blZeJjy3Ky8LHNuVl4eNd5WXhY6fysvCxS3lZ+NijvEn4
2Kd9JuHjgPIm4eOg8ibh47DyJuHjqPKy8fFX5WXj433lZds+p2CH8rLx0aW8
bHz0KC8HH6eVl4OPD5SXg4+PlJeDj4+Vl4OPT5SXi49PlZeLj8+Ul4uPz5WX
i48vlJeLj6+Ul4ePG8rLw8c3ysvDx7fKy8PHXeXl4eN75eXj4wfl5ePjgfLy
8fFIefn4eKK8fHz8rLwCfAw4afIK8DHopMkrcO9xf8hJk1eAj7STJq8AH8NP
mrxCfIxUXiE+RimvEB8ZyivEx1jlFdqjnAEnKG8yPrKUNxkfOcqbjI8C5U3G
R5HyJuOjRHlF+ChXXhE+qpRXhI9a5RXho0F5RfiwlFeMj4DyivERUV4xPuLK
K8ZHSnnF+HhReVPw8TvlTcHHS8qbgo9ZypuCjznKm4KPl5VXgo9XlFeCjwXK
K8HHq8orwcdryitx7/F4ifJK8dFkHp9Rio8VJ83fN6X4WGmenyrFx2q9Xik+
mk1ecyk+1un1S/GxweS1l+Jjk/qU4qPF5N0qxcdW9SvDx3b1K8PHDvUrw0er
+pXhY7f6leFjr/qV4WO/+pXho039yvBxSP3K8HFE/crwcUz9yvHRrn7l+Diu
fuX46FS/cnx0q185PnrVrxwfZ9SvHB9n1a8cH+fUr9ze6vi4oH7l+LiofhX4
uKR+Ffi4rD8/Ffi4qn4V+LimfhX4uK5+Ffi4qX4V+Lht8tor8HFH/SrwcU/9
KvBxX/28+PhR/bz4eKh+Xnw8Vj8vPp6qnxcfnk7Tz4uPgZ2mnxcfgztNPy8+
hnaafl58DOs0/bz4GGHyPJX4SDd5GZX4GG3yKirxMcbkpSrde3w5zuQtqMRH
pvpV4iNb/Srxkad+lfgoVL9KfBSrXyU+ytSvCh9e9avCR7X6VeGjTv2q8OFT
vyp8+NWvCh8h9avCR1T9qvCRUL8qfExVvyp8TFO/anzMUL9qfMxUv2p8zFa/
anzMVb9qfMxTv2p8zFe/anwsVL9qfCxSv2p8LFa/anwsU78afCxXvxp8vK5+
NfhYpX41+FijfjX4WKt+NfhYr341+NiofjX42Kx+NfjYon41+NimfrX4eFf9
avGxU/1q8bFL/WrxsUf9avGxT/1q8XFA/WrxcVD9avFxWP1q8XFU/Wrx8Z76
1eHjffWrw0eH+tXho0v96vDRo351+DitfnX4+ED96vDxkfrV4eNj9avDxyfq
V4ePT9WvHh+fqV89Pj5Xv3p8fKF+9fj4Sv3q7RuOjxvqV4+Pb9SvHh/fql89
Pu6qXz0+vle/enz8oH4N+Higfg34eKR+Dfh4on4N+PhZ/RrwMaDL9GvAx6Au
068BH0O6TL8GfKR1mX4N+BjeZfo14GOkyfM04mOUyctoxEeGyatoxMdYk5dq
xMcEk7egER9Z6teIjxz1a8RHgfo14qNI/RrxUaJ+je49PixXPx8+qtTPh49a
9fPho0H9fPiw1M+Hj4D6+fARUT8fPuLq58NHSv18+HhR/Xz4mK5+Fj5eUj8L
H7PUz8LHHPWz8PGy+ln4eEX9LHwsUD8LH6+qn4WP19TPwscS9bPw0aR+Nj5W
mMen2fhYqb42Plab5+fa+GjuMv++sfn8WGfybBsfG/R+bHxsMnmzbHy06P3Z
+Nhq8ppsfGzX+7XxscPktdj4aNX7t/Gx2+S12fjYa/Labf5/vt/k9dr4aNM+
Nj4OmbwrNj6OaC8bH8dMXp+Nj3bt58fHce3nx0en9vPjo1v7+fHRq/38+Dij
/fz4OKv9/Pg4p/38+Lig/fz4uKj9/Hx+XNJ+fnxc1n5+fFzVfn58XNN+fnxc
135+fNzUfn583NZ+fnzc0X5+fNzTfn583Nd+fv5/3qf9Avh4qP0C+His/QL4
eKr9AvjwnDL7BfAx8JTZL4CPwafMfgF8DD1l9gvgY9gps18AHyNMXlMAH+km
rzmAj9EmryWAjzEmrzWAj3Emry2Aj0yT1x7AR7bJ6w3gI8/knQ/go9DkXQng
o9jk3Qrgo8zk9QXw4TV5niA+qk1eWhAfdSYvI4gPn8nLDeLDf8r8/RPk8yOk
/YL4iGq/ID4S2i+Ij6naL4iPaSavKYiPGdoviI+Z2i+Ij9naL4iPudoviI95
Jq89iI/52i+Ij4XaL4iPRdoviI/F2i+Ij2Umry/I58dy7RfCx+vaL4SPVdov
hI812i+Ej7XaL4SP9dovhI+N2i+Ej83aL4SPLdovhI9t2i/E58c72i+Ej53a
L4SPXdovhI892i+Ej33aL4SPA9ovhI+D2i+Ej8PaL4SPo9ovhI/3tF8IH+9r
vzA+OrRfGB9d2i+Mjx7tF8bHae0X5t9XH2i/MD4+0n5hfHys/cL4+ET7hfHx
qfYL4+Mz7RfGx+faL4yPL7RfGB9fab8wPm5ovzA+vtF+YXx8q/3C+Lir/cL4
+F77hfHxg/YL4+OB9ovg45H2i+DjifaL4ONn7RfBx4Bus18EH4O6zX4RfAzp
NvtF8JHWbfaL4GN4t9kvgo+RJq8pgo9RJq85go8Mk9cSwcdYk9cawccEk9cW
wUeWyWuP4CPH5PVG8FFg8s5H8FFk8q5E8FFi8m5F8FFu8voi+Kg0eZ4oPmpN
XloUHw0mLyOKD8vk5UbxEdB+UXxEtF8UH3HtF8VHSvtF8fGi9oviY7r2i+Lj
Je0Xxccs7RfFxxztF8XHy9ovio9XtF8UHwu0XxQfr2q/KD5e035RfCzRflF8
NGm/KD5WaL8YPlZqvxg+Vmu/GD6atV8MH+u0XwwfG7RfDB+btF8MHy3aL4aP
rdovho/t2i+Gjx3aL4aPVu0Xw8du7RfDx17tF8PHfu0Xw0eb9ovh45D2i+Hj
iPaL4eOY9ovho137xfBxXPvF8dGp/eL46NZ+cXz0ar84Ps5ovzg+zmq/OD7O
ab84Pi5ovzg+Lmq/OD4uab84Pv6h/eL4uKr94vi4pv3i+Liu/eL4uKn94vi4
rf3i+Lij/eL4uKf94vi4r/3i+OjTfnF8PNR+CXw81n4JfDzVfgl8eHrMfgl8
DOwx+yXwMbjH7JfAx9Aes18CH8N6zH4JfIwweQsS+Eg3eU0JfIw2ec0JfIwx
eS0JfIwzea0JfGSavLYEPrJNXnsCH3kmrzeBj0KTdz6Bj2KTdyWBjzKTdyuB
D6/J60vgo9rkeZL4qDN5aUl8+ExeRhIffu2XxEdI+yXxEdV+SXwktF8SH1O1
XxIf07RfEh8ztF8SHzO1XxIfs7VfEh9ztV8SH/O0XxIf87VfEh8LtV8SH4u0
XxIfi7VfEh/LtF8SH8u1XxIff9J+Kdv9/irtl8LHGu2Xwsda7ZfCx3rtl8LH
Ru2Xwsdm7ZfCxxbtl8LHNu2Xwsc72i+Fj53aL4WPXdovhY892i+Fj33aL4WP
A9ovhY+D2i+Fj8M9vf8DpwFw2Q==
           "]]}, 
         {RGBColor[0.880722, 0.611041, 0.142051], PointSize[
          0.006944444444444445], AbsoluteThickness[1.6], 
          LineBox[{{0.9464308259088959, -4773.179999999999}, {
           1., -7325.377043011457}, {1.1, -6364.0952642008115`}, {
           1.144176656239496, -4773.179999999999}}], 
          LineBox[{{1.560712505966938, -4773.179999999999}, {
           1.6, -5888.08218312095}, {
           1.7000000000000002`, -6875.281973251477}, {
           1.8, -7072.271611710171}, {
           1.9000000000000001`, -5371.68895330659}, {
           1.9719616510711102`, -4773.179999999999}}], LineBox[CompressedData["
1:eJw9mXc8lf3/x900dBsRoRLZezsc+3Wcw5mVlDRokpYmTUnjTkTRUJoqaUkD
ZUSnpaF1hztpGEl72s2f+/f9vO/rjzzOo8/1/Dw/z+tyXcfjYzRtfnCkooKC
gm3PP//+jAvVVHQzUICwcYW2/ZFz8rt3/j0UEF+V7v7w9QX57sx/jz8QfNjw
7bBfhfKoGf8eitic4rXRq/q83NXl30MJHRq/Tt9WPy9X+P+jF6KPXm8cIrnA
eL1Q6ntd29D8POP1hvaPiMkfbEoZrw/8Mo0mNF4rY7y+EM6O+LNicDnjKePJ
wMT2I3w54ynDMKl83NVRcsbrh8lvBz7t317OeH+i/fuXslm8csZTwYtXM/vV
3ShjPFWYWl0eXKhUxniq6O/O+br6VynjqcE43KM9aX0J46mjeuSGK7vbihiv
P3a8aJpb8bmY8TRw4myzQDWJxmvgRvTHigk1JYyvgWcrtYrf/SiR///pLpr4
ytmoI2mn+TTRFiVtOPyuhPEGICsjeHrpCOozAI7RQSMTVYoZbwC8ReohXo/J
Rwvuokzp8yNFjKeFz1V6baedyU8b3rnXyoL2FzGeNhzCVbrHxRFPGx/LVv2V
GV/MeAMxSUdPtvMp8QZi782Iuo4LxNMBp9s9ZvX0YsbTwWNDPd3ThsTTQUhQ
xeSLHOLpQiy/0dUhLmY8XRgW/zL73EU8PdgdGf/i5mvi6cFC7mPh2UQ8PZwK
rlfWzCbeINS+GfUqUL2E8QYhrDz12bcE6jcYWVubVboy6HoMxseJ57s6ttD1
GIzHVupPyhuINwRNK3uf6/OU/IbgbdytsruPyU8fIza8uD77IPnp45dyjuf9
KvLTh9UbR5eUv4g3FIdWTd/V7wTxhiJgaqL6jmvEM8DOtMKKIduJZwDhkVwv
syXEM0DapMr5TUuIZ4hOvw4NQxviGWKV6UzbtTHEG4YUncbuzwuL2f01DDN2
929YcYb4wxC8IrT41lw6fxjOzT1UkBVN8w1DlYWl9s0DxDOCTkRIyts1xDPC
CgOfv6ZsI54Rgk6OchylTP2NoOEX4nVPn/oa4dlVpRGVXLoexljzx/3+au50
/xuj4Jvy2IkCuj7GsFiQHrSSQzxjNH6LL3YLJJ4xrh1rs5nvQjwT7CgfdbGq
Z/7/9THBCr+X467rEc8EE/Lq1xhrE88EWtnGBZs7ab0m+Duo3W7NV1qvKYy+
FxySVVNvUzyIvXU56iyt1xTH1x5/632W+pni5IO4cYuuEs8U201/rkh/Tjwz
mPcO9ta8RTwz2OxPETcdIJ4ZYge0GwYmEs8M6muPJsQtI54ZugcatHbtJZ45
vsxtWFKYTDxzmH2p8glNIp45GsxKXC8vJJ45hs7ctvDScuKZY/rtI8f804ln
gajMF8sqZxDPAsbdWbEzvIlngf66WmofbYlngVKLB5O3WBLPAnfMkkNk+sSz
RNrRxTFzlYlnCWm/xNxvX+h5Y4lb2yZ83FBPzxNLBC12PGb1uojxLLHW3XTs
4uYixrPCPIdFKUf+e75ZIa7zoOW8SuJZYeECY8PzVcSzwve8zcLXxcSzQqjK
s/gZ54hnjTcO5fJ1l4lnjaqxsVfibxPPGvlbp1auuU88ayT/Ov5i9EviWcMk
9mLfsnbi2aCtKfVE/RPi2eC+0VY0fieeDV6nnjSc8pp4Nnj/JnrQ5jbi2SC2
f4ZZjhL1s8Wk+njbe1XEs4X2p+wm+yfEs4XVXM6WCY+IZ4tbj884av7Hs0WL
5i6euTLx7DBq2/fO7i7i2WFgjurw+v+uhx2qFaa/6POZeHaQxpV89PtKPDu0
a/PHc77Qeu2x45+a69rfiWePguzktuSfxLOHfdmpo7WD6X6xR+mHK2csdeh+
sUdc5MzNnGHk54CkrFen5+nT/eIAJ6G4wtWC7j8HxKoXHHK3JJ4D9v3tpbZj
NPEc8NJgod/cqcRzhPLvxkF1YWz8IEesTtI5mB9Kzy9HbNKfrhY4jp0vc4Rp
2FfDyxE0X8//c/Xn3F7AeAmOWN+yYFlSLM3viEbVDoM38YyX3zM+sYhvuYJ8
HBEY0Ka8PIbxWhyRdc7i+7xY8nPCyzENoTNXkp8Tnl1ptXm9nPycsFeWfEwW
Q35OKMyepqiygvyccKJvuKJmHPk54e4C7y1LFpCfEw6VVSSvXUl+Togs3n8i
eA35OSEpuOaoZAP5OeHFAMG0qAzyc8bKOLeFCYfJzxlpmpdsS06RnzPEO1vz
hp4iP2eoxleV7i8gP2eMHKkUOPYy+TnjkPyk39Db5NczXjPm7eRL5OeMBe8t
dhaXkJ8zNv8dfe1nKfk5IzayS7b5JPm5IE474ZBWIfm5IMLyn2kGRXT/uKCZ
uzQ+r4D8XDDrzNKKfOLNcIH7qq2TD5aSnwvC+zp3Si+Rnwv48kuSw5WMl++C
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
5uDRcvn/AYebBRo=
           "]], 
          LineBox[{{0.06803168781795584, -4773.179999999999}, {
           0.0747895049220504, -9203.900000000001}}], 
          LineBox[{{0.16844902390431823`, -9203.900000000001}, {
           0.18679717546060257`, -4773.179999999999}}]}}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0, -4551.644},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Speed / [foot/sec]\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Speed Truth and Estimate\.7f vs. Time\"", FontFamily -> 
            "Candara", FontSize -> 12, StripOnInput -> False], 
           TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0, 57.5}, {-9203.900000000001, -4773.179999999999}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.05], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}]}]}
   },
   AutoDelete->False,
   GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}],
  "Grid"]], "Output",
 CellChangeTimes->{
  3.6699227994811573`*^9, 3.669922853786981*^9, 3.669922888960762*^9, {
   3.669922977116374*^9, 3.6699229885365047`*^9}, 3.669923281031231*^9, 
   3.669923327748118*^9, 3.669923875743981*^9, 3.669923919912553*^9, 
   3.669923955267486*^9, {3.669924006575589*^9, 3.669924036568014*^9}, 
   3.66992434377634*^9, 3.669924386397457*^9, 3.6699246697373962`*^9, {
   3.6699247064136677`*^9, 3.669924736356739*^9}, {3.6699248123832197`*^9, 
   3.669924820693666*^9}, 3.669924854365623*^9, 3.6699251085376177`*^9, 
   3.669925364221159*^9, 3.669925681898447*^9, 3.669925768191081*^9, 
   3.669925902589381*^9, 3.669926001315505*^9, {3.669926085091527*^9, 
   3.66992609986534*^9}, 3.669926164721451*^9, {3.669926205358274*^9, 
   3.6699262163847027`*^9}, 3.669926448891405*^9, 3.6699265697910833`*^9, {
   3.669926637235973*^9, 3.669926659729146*^9}, {3.669926706840378*^9, 
   3.6699267587051563`*^9}, 3.669926835635337*^9, 3.6699268783318*^9, {
   3.669926950163229*^9, 3.669927075141576*^9}, {3.6699271102637997`*^9, 
   3.669927126965555*^9}, 3.669929261774242*^9, 3.669929293622839*^9, 
   3.669929347844597*^9}]
}, Open  ]],

Cell["\<\
Doesn\[CloseCurlyQuote]t work well for a two-state system, but we can fix \
this (later).\
\>", "Text",
 CellChangeTimes->{{3.6699275304863443`*^9, 3.669927545284595*^9}, {
  3.66992793649763*^9, 3.669927938333159*^9}}],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"With", "[", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{"x0", "=", "400000"}], ",", 
     RowBox[{"v0", "=", 
      RowBox[{"-", "6000"}]}], ",", 
     RowBox[{"a0", "=", 
      RowBox[{"-", "32.2"}]}], ",", 
     RowBox[{"\[CapitalZeta]", "=", 
      RowBox[{"col", "[", 
       RowBox[{"{", 
        SuperscriptBox["1000.0", "2"], "}"}], "]"}]}]}], "}"}], ",", 
   "\[IndentingNewLine]", 
   RowBox[{"experiment", "[", "\[IndentingNewLine]", 
    RowBox[{"0", ",", "30.", ",", "0.1", ",", "\[IndentingNewLine]", 
     RowBox[{"col", "[", 
      RowBox[{"{", 
       RowBox[{"x0", ",", "v0"}], "}"}], "]"}], ",", "\[IndentingNewLine]", 
     RowBox[{"1000000", 
      RowBox[{"id", "[", "2", "]"}]}], ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"dt", ",", "t"}], "}"}], "\[Function]", "\[IndentingNewLine]", 
      
      RowBox[{"{", 
       RowBox[{
        RowBox[{"1000000", "*", 
         RowBox[{"(", GridBox[{
            {
             FractionBox[
              SuperscriptBox["dt", "3"], "3"], 
             FractionBox[
              SuperscriptBox["dt", "2"], "2"]},
            {
             FractionBox[
              SuperscriptBox["dt", "2"], "2"], "dt"}
           }], ")"}]}], ",", "\[IndentingNewLine]", 
        RowBox[{"(", GridBox[{
           {"1", "dt"},
           {"0", "1"}
          }], ")"}], ",", "\[IndentingNewLine]", 
        RowBox[{"(", GridBox[{
           {"dt"},
           {"1"}
          }], ")"}], ",", "\[IndentingNewLine]", 
        RowBox[{"col", "[", 
         RowBox[{"{", "g", "}"}], "]"}], ",", "\[IndentingNewLine]", 
        RowBox[{"(", GridBox[{
           {"1", "0"}
          }], ")"}], ",", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"col", "[", 
          RowBox[{"{", 
           RowBox[{"x0", "+", 
            RowBox[{"v0", " ", "t"}]}], "}"}], "]"}], "+", 
         RowBox[{"gen", "[", "\[CapitalZeta]", "]"}]}]}], "}"}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"t", "\[Function]", 
      RowBox[{"x0", "+", 
       RowBox[{"v0", " ", "t"}], "+", 
       RowBox[{"a0", " ", 
        RowBox[{
         SuperscriptBox["t", "2"], "/", "2"}]}]}]}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"t", "\[Function]", 
      RowBox[{"v0", "+", 
       RowBox[{"a0", " ", "t"}]}]}]}], "]"}]}], "]"}]], "Input",
 CellChangeTimes->{{3.6699240792959423`*^9, 3.6699240828557367`*^9}, {
   3.6699241417974777`*^9, 3.669924176659999*^9}, 3.669924242073472*^9, {
   3.669924400254054*^9, 3.669924402428566*^9}, {3.669924467961753*^9, 
   3.669924626234639*^9}, {3.669925800200728*^9, 3.669925810731179*^9}}],

Cell[BoxData[
 TagBox[GridBox[{
    {
     RowBox[{
      GraphicsBox[{{}, {{}, {}, 
         {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
          0.011111111111111112`], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1: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
US76LEih/0HrvQ98bhak+rn/YpMxuejj9TKZ6YthSvQxr1kZmuEKXxKvd7hw
IPVYDvp43inI8tzVXPQxz0vz/foRvkbmyUHdtKfAJyyJdrds3PdYiT5L3p9a
VjP1LO6PJalWGMjVunno4/mZbxWrvc+gj/drB+i0e3AWfZakcJ+dqxp2Fn28
f8f1wrWleT98OTyv0Bnb6AJfFe9fu3hDatc89DH3SytN7ZSHPiui0r8ejGs6
gz4rUmW80vMNVqKPuZ+fkVMm3tfNioTrzeLmdVIfc8vTSUHdpT7mhJC/uqyX
+thfMLHj9HipjzngjZHhuHz0WZEiuXCf0jsffTyPd00JX38OfdYkRgXLco7n
o8+a79vNlWY9zqOPeWhL88+R8LlZk6p+nIGD93n0MVd/izDqdQF9vF7f5Y/R
2y+gz5ooqn3zqLgC9PH6z0tnlbgUoo/P76L3U/n6AvRZk0dn2WH92YXos6GG
N90qStYVos+GxJWj055+KUSfDRmke9dELCtCnw2pFm5XbmsqRB+vd6wLD2mG
T8HzUwvaVndXoc+G1C0jlGtnqdDH8xtbb7f8Bl8Vc8hvPU9cgq/Rhig4r+P9
hUXosyXF7WLLL3OK0GdLHvuO6VmnwmdpS2adE4w3LVGhj1nzbvkDDxX6bKnh
+cohemHwKWwp69Eum+BdUp8tpaQb+q27U4Q+Pi+/dfitWqnPlnxNXdLv+8PX
yH6Hh5PsElTosyOxeHq2kU8x+uyoemyJ463RxeizI7W6X53PpmL02ZHBxasH
1n8qRp8d+X47o597oQR9dpQlDtQXvS9Bnx012Gz2L/YpQR/zoa91M+JK0Mfn
O8Vo7pp5EX18fs89PT8cuog+e5IvKy/Y0P0i+uzJ4/NeH12Ni+izp2oqyBqv
uog+5l8+Fn4OLkWfPdG0zh6lRqXosyftSd5DTgSWos+eUk7pBXY9WYo+e2ow
qA3aV16KPnvy7d3l86bkMvTx+hc6K+MMy9DnQNrXno7bubQMfQ7UoJrjOdO+
DH0OFNdi9GdaXhn6HKh6iXHvhLnl6OO52qVOK7scfQ6k9pr8eyebS+hzILLx
elKSVY4+B5KHB0z/lFuOPgdSPU8p7vkYvkbuOe3vKr8Hn3Aks3+y21JaOfoc
SX3NKSrd6xL6HCllUsDdaXGX0OdIpFVweOTRS+hzJN/qJr21Oy+hz5EaqnWL
XvWS+hxJUTXNc8c4+HKYqx0iqyLgq3Ik7bHzZOoW+BqZFc0TYqvgE05UHfrh
1LUeFehzIsXR9730lBXoc6KUhZYrVvtdRp8TqZ/0zGj38jL6eP21Q2lpiVfQ
50RxMUEzVk+pRB/Plaf6yxuvoM+J5EvbrDlbfwV9TmRWMe+q3L0SfU6UdbRU
duFeJfqcyaAoONvo5RX0OZPaojir+St8ls7kMeZbU59U+NycSfvnJ6Fzv1xB
H/P5B51nRsCncCZf+buBWVulPmdKiVxrvv5qJfrYf17LYVJTJfqcSV6VXrsu
vAp9zkRGGyY9M61C3xiSNy2qX51dqfovEtQ14Q==
           "]]}}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0, 0},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Position Residual / foot\"", FontFamily -> "Candara", FontSize -> 
            12, StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Position Residuals vs. Time\"", FontFamily -> "Candara", 
            FontSize -> 12, StripOnInput -> False], TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0, 30.}, {-14598.897114610445`, 165.9576060026884}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.05], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}], " ", 
      GraphicsBox[{{}, {{}, {}, 
         {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
          0.011111111111111112`], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1: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
Xm8KSNtdNRhNPNkHeB5qpV3q/wCMGC+h
           "]]}, 
         {RGBColor[0.880722, 0.611041, 0.142051], PointSize[
          0.011111111111111112`], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1: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
           "]]}}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0, 195778.12529753562`},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Position / foot\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Position Truth and Estimate vs. Time\"", FontFamily -> 
            "Candara", FontSize -> 12, StripOnInput -> False], 
           TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0, 30.}, {205510., 400147.4940492877}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.05], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}]}], 
     RowBox[{
      GraphicsBox[{{}, {{}, {}, 
         {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
          0.011111111111111112`], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1: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
8SXk54G6laHRKxJKyI8br2YfGs8i3lUPfOl/ojqmaCv5cfNnH9fvXllMflKY
zFk/03o47WdTKWwc4/GfHfz3mxRzNvu6Ox2n75V5UlQaLAztE1eSnxShJwxO
/+cLPn9SHDlblBvUU0V+Ukhex898vJ8/Xzj+jHP2IR6F5CdF6zsLgcHBQvKT
QhM91DL6PEt+nmDPzFj9tDaR/DyheZfj5bWR719PdA2T97XNLCY/T+gmq1fW
lRMv0hMJ16S1s6/nkB/HUyn6rYxCyc8TjNv79zeZSPLj4iW7bkj9EsnPE6JW
30bJD9RvnZ5QLXmkSexLIj8vKI9eyNhwJpP8vNC6x7ct7yit18ULp+b+clE/
fgv5eaEwJiS0sa6K/Lyw4oNl3+2sb8jPC1n1NX79njXk54Wy3S3ZOeFV5OcF
sUJ3SKAHf75wPPmdxYPH+fp6obLwN5f3wmry84aopmTZpkD6fjL15vrjxm7/
UP7/A95gHScMev4ZSX7eELJ3VzhvSyU/bzStzW4ePFtAft5A0H2B8Fvaz9u9
oWxMv7fMl+8Pb8z5+tKUS+D7g4vrqg17Z5aSnzfk/j/WL4jj6/s5cEazZnj/
StX/ABotGRU=
           "]]}}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0, 0},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Speed Residual / [foot/sec]\"", FontFamily -> "Candara", 
            FontSize -> 12, StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Speed Residuals vs. Time\"", FontFamily -> "Candara", FontSize -> 
            12, StripOnInput -> False], TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0, 30.}, {-1742.9973434755511`, 1082.6265432773416`}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.05], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}], " ", 
      GraphicsBox[{{}, {{}, {}, 
         {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
          0.008333333333333333], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1:eJxV1HlM03ccxvHKYDKGAxHEIlc5W87ed/n9NGZshrjgzMwicXFzmIVFMTrj
FYPxiFc84jQ6o9EwDcZI3IhSaDANxmM6jdMYnUe84jHRyJyoYzrX7/f5Pn+s
f2jE9t2nfHhh+nLWxK/jDAZD5RCDQfz9v0dzZ/SH7eLRpeWt6Kv3zuG/f9FC
SZdbm+d1Rp0O8bioNWzoHWhbwP//XVuU3l57e3FnFKFbmvi6sYXPv6uFc2PB
ZZ3RGY3icV+L1UJrVvL1f2gvzLHi6s7o2V/F47EmcoPr2OvXZGYjn/+nNjEW
bNrM/l/abDFwK1//XIvNM1/fzvd7ocVyC9N3svdKE8+r283e39rj+liwlb1/
NPFxe/ax91qziA+8n71/tY9iuaqD7Bl0+eVDfL5Bj33a8K4O9ofoPzbHBh7h
6+P0YwPmhpRuvt87eizXXtvDXrwu/myJshevi3uEj7GXoIuP23+CvXd1Mc98
mr2hurjHtLPsJerydefZS9TFPX67yN57urhH0mX2knRxj7FX2XtfF/dYeIO9
ZF18teMWe8m6uEffXfaG6eIehQ/Y+0AX8xoesZeii9z3T9hL1WWnn89P1cU9
Ep6zn6qLe4Redkblyx3DdXGP7wb5fsN1cY/2N+ylyXvcN4RVL03eIzc+rHpp
8h6Th4bVnhHyHhuSwqo3Qt7j5LCw6qXLexiGs5cu7+FNZy9d3qM5k70MeY+2
LPYy5D1u57A3Ut7DaGJvpLxHfRF7I+U91pSylynv0VvGXqa8x2Ale6PkPRw2
9kbJezQ52Rsl79HqYc8o73Hdz55R3iO9hr0seY+6MexlyXusGMdelrxHTy17
o+U9BsazN1reo2oCe9nwUc9eNnxMYi8bPiazlwMfU9jLgY+p7OXCxzT2cuFj
Onu58DGDvTz4aGIvDz5mspcPH7PD6ucrHz7msp8PH/P5+nz4WMT3y4ePJeyZ
4GMpeyb4WM6eCT5WsWeCj7XsmeBjPXsF8LGJvQL42MJeAXxsY68APnawVwAf
u9grhI89/P4Uwsde9grho429Qvg4wF4hfLSzVwQfP7NXBB+H2SvSa8TAMHtF
8BFhrwg+jrJXDB+97BXDx3H2iuHjFHvF8HGGvWL4OMdeCXxcYK8EPi6xVwIf
V9grgY9r7JXAx032SuHjDnul8HGPvVL4eMheKXz0sVcKH0/ZM8PHM/bM8DHA
nhk+XrFnho/X7Jnh4y17FviI61I9C3wkdKmeRd6jP7FL9Szwkdylehb4SOlS
vTL4SGOvDD4y2CuDDyN7ZfCRzV6ZniG+gXnslcNHAXvl8FHMXjl8WNgrh48K
9srho5q9Cviws1cBH272KuDDx14FfATZq4APjb1K+BjLXiV8fMheJXx8zF4l
fNSxVwkfn7BXBR+fslcFH5+xVwUfn7NXBR8N7FXBxxfsVcPHV+xVw0cje9Xw
8Q171fDxLXvV8h6Ds9izwscc9XyjFT7mdanfN1b4WKBeX2eFj8V8Pyt8tKhe
ixU+lvH9rfCxUvU6rPCxmnus8LFO9e5b4WMj99ngYzP32eBjK/fZ4GM799ng
Yyf32eBjN/fZ4KOV+2zwsY/7bPCxn/ts8HGQ+2zwcYj77PDRwX12+DjCfXb4
6OY+O3z0cJ8dPqLcZ4ePY9xnh48T3GeHj9PcZ9c3Ch9nuc8OH+e5zwEfF7nP
AR+X+fPjgI+r3OeAjxvc54CPW9zngI+73OeAjweq1+GAj0fc54CPJ9zngI9+
7nPCx3Puc8LHS+5zwscg9znh4w33OeHD0K32OeEjvlvtc8LH0G61zwkfSd1q
nxM+hnWrfU74SFU9gws+0lXP6IKPTNVzuOAjS/XqXPIe13NUr9EFHybuc8FH
Efe54KOU+1zwUcZ9Lvio5D4XfNi4zw0fTu5zw4eH+9zw4ec+N3zUcJ8bPsZw
nxs+xnGfGz5quc8NH+O5zw0fE7jPDR/13OeBj0nc54GPydzngY8p3OeBj6nc
54GPadzngY/p3OeBjxnc54GPJu7zwMdM7vPAx2zu88LHXO7zwsd87vPCxyLu
88LHEu7zwsdS7vPCx3Lu88LHKu7zwsda7vPCx3ru88LHJu7zwccW7vPBxzbu
88HHDu7zwccu7vPBxx7u88HHXu7zwUcb9/ng4wD3+eCjnft88PET9/nh4zD3
+eEjzH1++Ihwnx8+jnKfHz56uc8PH8e5zw8fp7jPDx9nuM8PH+e4zw8fF7gv
AB+XuC8AH1e4LwAf17gvAB83uS+g3xE+7nBfAD7ucV8APh5yXwA++rgvAB9P
uS8AH8+4LwgfA9wXhI9X3BeEj9fcF4SPt9wXhI+4iNoXhI+EiNoXhI/EiNoX
hI/kiNoXhI+UiNoXhI801TOE4CND9Ywh+DCqniMEH9mqVxeCjzzVawzBRwH3
heCjmPtC8GHhvhB8VHBfCD6quS8k73HSzn018OGORP8DvUSHKg==
           "]]}, 
         {RGBColor[0.880722, 0.611041, 0.142051], PointSize[
          0.008333333333333333], AbsoluteThickness[1.6], 
          LineBox[CompressedData["
1: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
r5D/FzgLqWI=
           "]]}}, {}},
       AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
       Axes->{True, True},
       AxesLabel->{None, None},
       AxesOrigin->{0, -5059.334224994658},
       DisplayFunction->Identity,
       Frame->{{True, True}, {True, True}},
       FrameLabel->{{
          FormBox[
           StyleBox[
           "\"Speed / [foot/sec]\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox["\"\"", TraditionalForm]}, {
          FormBox[
           StyleBox[
           "\"time / sec\"", FontFamily -> "Candara", FontSize -> 12, 
            StripOnInput -> False], TraditionalForm], 
          FormBox[
           StyleBox[
           "\"Speed Truth and Estimate\.7f vs. Time\"", FontFamily -> 
            "Candara", FontSize -> 12, StripOnInput -> False], 
           TraditionalForm]}},
       FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
       GridLines->{Automatic, Automatic},
       GridLinesStyle->Directive[
         GrayLevel[0.5, 0.4]],
       ImagePadding->All,
       ImageSize->Medium,
       Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& ), "CopiedValueFunction" -> ({
             (Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
              Part[#, 1]], 
             (Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
              Part[#, 2]]}& )}},
       PlotRange->{{0, 30.}, {-7482.071287120274, -5174.702656524449}},
       PlotRangeClipping->True,
       PlotRangePadding->{{
          Scaled[0.02], 
          Scaled[0.02]}, {
          Scaled[0.05], 
          Scaled[0.05]}},
       Ticks->{Automatic, Automatic}]}]}
   },
   AutoDelete->False,
   GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}],
  "Grid"]], "Output",
 CellChangeTimes->{
  3.66992419077186*^9, 3.669924516638192*^9, {3.669924555238607*^9, 
   3.669924627803458*^9}, 3.669925378448237*^9, 3.669925812162897*^9}]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1075, 855},
WindowMargins->{{0, Automatic}, {Automatic, 0}},
FrontEndVersion->"10.4 for Mac OS X x86 (32-bit, 64-bit Kernel) (February 25, \
2016)",
StyleDefinitions->Notebook[{
   Cell[
    StyleData[
    StyleDefinitions -> 
     FrontEnd`FileName[{"Report"}, "StandardReport.nb", CharacterEncoding -> 
       "UTF-8"]]], 
   Cell[
    StyleData["Title"], FontFamily -> "Candara", FontSize -> 44, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Subtitle"], FontFamily -> "Candara", FontSize -> 24, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Subsubtitle"], FontFamily -> "Candara", FontSize -> 16, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Author"], FontFamily -> "Candara", FontSize -> 14, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Department"], FontFamily -> "Candara", FontSize -> 11, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Date"], FontFamily -> "Candara", FontSize -> 11, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Chapter"], FontFamily -> "Candara", FontSize -> 34, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Subchapter"], FontFamily -> "Candara", FontSize -> 28, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Section"], FontFamily -> "Candara", FontSize -> 28, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Subsection"], FontFamily -> "Candara", FontSize -> 20, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Subsubsection"], FontFamily -> "Candara", FontSize -> 19, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Text"], FontFamily -> "Candara", FontSize -> 15, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Item"], FontFamily -> "Candara", FontSize -> 15, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["ItemParagraph"], FontFamily -> "Candara", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Subitem"], FontFamily -> "Candara", FontSize -> 15, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["SubitemParagraph"], FontFamily -> "Candara", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Subsubitem"], FontFamily -> "Candara", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["SubsubitemParagraph"], FontFamily -> "Candara", FontSize -> 15,
     FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["ItemNumbered"], FontFamily -> "Candara", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["SubitemNumbered"], FontFamily -> "Candara", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["SubsubitemNumbered"], FontFamily -> "Candara", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["DisplayFormula"], FontFamily -> "Candara", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["DisplayFormulaNumbered"], FontFamily -> "Candara", FontSize -> 
    15, FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["InlineFormula"], FontFamily -> "Palatino", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Code"], FontFamily -> "Inconsolata", FontSize -> 15, 
    FontWeight -> "Bold", FontSlant -> "Plain"], 
   Cell[
    StyleData["Program"], FontFamily -> "Courier", FontSize -> 15, FontWeight -> 
    "Plain", FontSlant -> "Plain"], 
   Cell[
    StyleData["Input"], FontFamily -> "Inconsolata", FontSize -> 15, 
    FontWeight -> "Bold", FontSlant -> "Plain"], 
   Cell[
    StyleData["Output"], FontFamily -> "Inconsolata", FontSize -> 15, 
    FontWeight -> "Plain", FontSlant -> "Plain"]}, WindowSize -> {718, 855}, 
  WindowMargins -> {{Automatic, 0}, {Automatic, 0}}, Visible -> False, 
  FrontEndVersion -> 
  "10.4 for Mac OS X x86 (32-bit, 64-bit Kernel) (February 25, 2016)", 
  StyleDefinitions -> "Default.nb"]
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1486, 35, 127, 2, 92, "Title"],
Cell[1616, 39, 99, 1, 43, "Author"],
Cell[1718, 42, 97, 1, 44, "Date"],
Cell[CellGroupData[{
Cell[1840, 47, 98, 1, 68, "Chapter"],
Cell[1941, 50, 301, 7, 82, "Input"],
Cell[2245, 59, 3372, 92, 395, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[5654, 156, 232, 3, 68, "Chapter"],
Cell[5889, 161, 215, 4, 34, "Text"],
Cell[6107, 167, 1987, 49, 154, "Input"],
Cell[CellGroupData[{
Cell[8119, 220, 1733, 48, 155, "Input"],
Cell[9855, 270, 4159, 112, 59, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[14051, 387, 1232, 26, 155, "Input"],
Cell[15286, 415, 12296, 308, 494, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[27619, 728, 138, 1, 63, "Subchapter"],
Cell[CellGroupData[{
Cell[27782, 733, 936, 20, 40, "Input"],
Cell[28721, 755, 2005, 52, 97, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[30763, 812, 939, 20, 40, "Input"],
Cell[31705, 834, 2000, 52, 99, "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[33766, 893, 106, 1, 68, "Chapter"],
Cell[CellGroupData[{
Cell[33897, 898, 108, 1, 63, "Subchapter"],
Cell[34008, 901, 437, 8, 57, "Text"],
Cell[CellGroupData[{
Cell[34470, 913, 1197, 35, 210, "Input"],
Cell[35670, 950, 1328, 39, 128, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[37047, 995, 97, 1, 63, "Subchapter"],
Cell[CellGroupData[{
Cell[37169, 1000, 406, 7, 28, "Item"],
Cell[37578, 1009, 481, 10, 28, "Item"],
Cell[38062, 1021, 553, 12, 28, "Item"],
Cell[38618, 1035, 498, 10, 28, "Item"],
Cell[39119, 1047, 537, 10, 28, "Item"],
Cell[39659, 1059, 494, 6, 28, "Item"],
Cell[40156, 1067, 563, 16, 28, "Item"]
}, Open  ]],
Cell[CellGroupData[{
Cell[40756, 1088, 6153, 157, 584, "Input"],
Cell[46912, 1247, 2123, 54, 97, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[49084, 1307, 112, 1, 63, "Subchapter"],
Cell[49199, 1310, 1123, 28, 57, "Text"],
Cell[50325, 1340, 153, 2, 34, "Text"],
Cell[50481, 1344, 4845, 123, 415, "Input"],
Cell[55329, 1469, 13823, 288, 1287, "Input"],
Cell[CellGroupData[{
Cell[69177, 1761, 6170, 135, 678, "Input"],
Cell[75350, 1898, 1853, 29, 75, "Message"],
Cell[77206, 1929, 2230, 35, 75, "Message"],
Cell[79439, 1966, 1856, 29, 75, "Message"],
Cell[81298, 1997, 2230, 35, 75, "Message"],
Cell[83531, 2034, 44984, 805, 573, "Output"]
}, Open  ]],
Cell[128530, 2842, 229, 5, 34, "Text"],
Cell[CellGroupData[{
Cell[128784, 2851, 2650, 72, 441, "Input"],
Cell[131437, 2925, 31408, 582, 531, "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}
]
*)

(* NotebookSignature bwDZGD1Q9FJ5#Ag530tzM#iU *)
