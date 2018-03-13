(* ::Package:: *)

(* Modules for topology optimization of 2D rectangular structures*)
(* by Luis Bahamonde Jacome *)
(* Copyright 2018  *)


BeginPackage["TopologyOptimzation2D`"]

TOP2D::usage = "TOP2D provides functions to support topology optimization for minimum compliance of 2D rectangular structures"

BisectionRoot::usage = "bisection method for solving roots of equations"

LocalUpdateRule::usage = "LocalUpdateRule[designData, simData, utilizationMetric] Setoodeh et al update rule"

DesignUpdate::usage = "update the current design to a new design"

InitialDesignVector::usage = "InitialDesignVector[modelData, initialValue]"

OptimalityCriteria::usage = "OptimalityCriteria[modelData, designData] Fixed-point iteration of the local update rule"



Begin["`Private`"]
Needs["SquareFiniteElement`", "fem2d.wl"]


(*  low level interface *)

LocalUpdateRule[designVector_List, penal_, Ymin_, strainEnergyDensities_List, utilizationMetric_] := Module[{newDensities},
newDensities = Surd[strainEnergyDensities / (utilizationMetric/penal), penal + 1]; (* nth root *)
newDensities = designVector * newDensities;
newDensities[[Position[newDensities,_?(#<Ymin&)] //Flatten]] = Ymin;
newDensities[[Position[newDensities,_?(#>1.&)] //Flatten]] = 1.;
newDensities
]

(* high level interface *)
LocalUpdateRule[designData_Association, simData_Association, utilizationMetric_] := 
LocalUpdateRule[designData["designVector"], designData["penal"], designData["voidModulus"], simData["strainEnergyDensity"], utilizationMetric]
                                                                                                     

(* it is implicit that the current design is in the designData data structure*)
DesignUpdate[modelData_Association, designData_Association, utilizationMetric_] := Module[{simData, designVector},
simData = SimulateFEModel[modelData, designData];
designVector = LocalUpdateRule[designData, simData, utilizationMetric];
{designVector, simData}
]

(* the user explicitly gives the current design*)
ExplicitDesignUpdate[modelData_Association, designData_Association, utilizationMetric_NumericQ, designVector_List] := Module[{currentDesignData},
currentDesignData = designData;
currentDesignData["designVector"] = designVector;
DesignUpdate[modelData, currentDesignData, utilizationMetric]
]

(* TODO: if statement to set the numElems or numNodes as an Option, by default elems *)
InitialDesignVector[modelData_, initialValue_:1.]:= ConstantArray[initialValue, modelData["numElems"]] 


(*BUG: when OptimalityCriteria calls this function it always converges to all zeros!!!!*)
FixedPointSolver[residual_, initialGuess_, tol_:1*^-3]:= FixedPoint[residual, initialGuess, SameTest -> (Max[Abs[#1-#2]]<tol&)]



(* volume fraction constraint *)

(* low level interface *)
VolumeFractionResidual[designVector_List, volumeFraction_] := Mean[designVector] - volumeFraction

(* high level interface *)
VolumeFractionResidual[designData_Association] := VolumeFractionResidual[designData["designVector"], designData["volumeFraction"]]

LagrangeMultiplierResidual[simData_Association, designData_Association, mu_] := Module[{designVector},
designVector = LocalUpdateRule[designData, simData, mu];
VolumeFractionResidual[designVector, designData["volumeFraction"]]
]

Options[BisectionRoot] = {"LowerBound" -> 1*^-9, "UpperBound" -> 1*^9}
BisectionRoot[residual_,  tol_:1*^-3, OptionsPattern[]] :=Module[{iterNum, midPoint, lower, upper},
iterNum = 0;
lower = OptionValue["LowerBound"];
upper = OptionValue["UpperBound"];
While[(upper - lower)/(lower + upper) > tol, 
iterNum ++;
midPoint = 0.5 * (upper + lower);
If[residual[midPoint] > 0., lower = midPoint, upper = midPoint]];
   {midPoint, iterNum}
]

DesignUpdate[modelData_Association, designData_Association] := Module[{simData, mu, designVector},
simData = SimulateFEModel[modelData, designData];
mu = BisectionRoot[LagrangeMultiplierResidual[simData, designData, #]&][[1]];
designVector = LocalUpdateRule[designData, simData, mu];
{designVector, simData, mu}
]


ExplicitDesignUpdate[modelData_Association, designData_Association, designVector_List] := Module[{currentDesignData},
currentDesignData = designData;
currentDesignData["designVector"] = designVector;
DesignUpdate[modelData, currentDesignData]
]

Options[OptimalityCriteria] = {"DensityFilter" -> Identity}
OptimalityCriteria[modelData_Association, designData_Association, OptionsPattern[]] := Module[{residual}, 
residual = OptionValue["DensityFilter"]@ExplicitDesignUpdate[modelData, designData, #][[1]]&;
(* FixedPointSolver[residual, designData["designVector"]] *)
FixedPoint[residual, designData["designVector"], SameTest -> (Max[Abs[#1-#2]]<1*^-2&)]
]



(* TODO: ParametricDesignStudy 
PrepareDesignRun[paramRule_, inputData_, designData_] := 
If[StringMatchQ[Keys[paramRule], #]& /@ ]
DesignRun[paramRule_, inputData_, designData_]:= Module[{actualInputData, actualDesignData, modelData},
{actualInputData, actualDesignData} = PrepareDesignRun[paramRule, inputData, designData];
modelData = ModelMBBBeam[actualInputData, actualDesignData];
OptimalityCriteria[modelData, actualDesignData]
]
ParametricDesignStudy[paramName_String, {start_, end_, step_}, inputData_, designData_]:=
Table[DesignRun[paramName \[Rule] thisValue, inputData, designData], {thisValue, start, end, step}]
*)


End[]
EndPackage[] 
