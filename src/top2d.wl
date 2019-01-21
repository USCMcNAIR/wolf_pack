(* ::Package:: *)

(* Modules for topology optimization of 2D rectangular structures*)
(* by Luis Bahamonde Jacome *)
(* Copyright 2018  *)


BeginPackage["TopologyOptimzation2D`"]

TOP2D::usage = "TOP2D provides functions to support topology optimization for minimum compliance of 2D rectangular structures:
UpdateKey
OptimalityCriteria
InitialDesignVector"

UpdateKey::usage = "UpdateKey[data, key, value] updates an association's key string with the given value"

BisectionRoot::usage = "bisection method for solving roots of equations"

LocalUpdateRule::usage = "LocalUpdateRule[designData, simData, utilizationMetric] updates the designVector in designData using Setoodeh et al. rule
LocalUpdateRule[designData, simData, utilizationMetric, \!\(\*
StyleBox[\"nodalField\",\nFontWeight->\"Plain\"]\)] updates the nodal designVector if nodalField is True
LocalUpdateRule[designVector, penal, Ymin, strainEnergyDensities, utilizationMetric] is the low-level interface"

DesignUpdate::usage = "update the current design to a new design"

InitialDesignVector::usage = "InitialDesignVector[modelData, initialValue]"

OptimalityCriteria::usage = "OptimalityCriteria[modelData, designData] Fixed-point iteration of the local update rule"



Begin["`Private`"]
Needs["SquareFiniteElement`", "fem2d.wl"]


UpdateKey[data_Association, key_?StringQ, value_]:= Module[{dataCopy},
dataCopy = data;
dataCopy[key] = value;
dataCopy
]


(* Nodes to Elements Mapping*)

NodalHomogenization[field_, penalization_, nodes_]:= Surd[Mean[Power[field, penalization][[nodes]]], penalization]
MapToNodes[field_List, modelData_Association]:= NodalHomogenization[field, modelData["designData", "penal"], #]& /@ modelData["nodalConnectivity"]
MapToNodes[field_?StringQ, simData_Association]:= MapToNodes[simData[field], simData["modelData"]]

ElemHomogenization[field_, penalization_, elems_]:= Surd[1/Mean[1/Power[field, penalization][[elems]]], penalization]
MapToElems[field_List, modelData_Association]:= ElemHomogenization[field, modelData["designData", "penal"], #]& /@ modelData["elemConnectivity"]
MapToElems[field_?StringQ, simData_Association]:= MapToElems[simData[field], simData["modelData"]]


(* Fully Utilized Design*)

LocalUpdateRule[designVector_List, penal_?NumericQ, Ymin_?NumericQ, strainEnergyDensities_List, utilizationMetric_?NumericQ] := Module[{newDensities},
newDensities = Surd[strainEnergyDensities / (utilizationMetric/penal), penal + 1]; (* nth root *)
newDensities = designVector * newDensities;
newDensities[[Position[newDensities,_?(#<Ymin&)] //Flatten]] = Ymin;
newDensities[[Position[newDensities,_?(#>1.&)] //Flatten]] = 1.;
newDensities
]

LocalUpdateRule[designData_Association, simData_Association, utilizationMetric_?NumericQ, nodalField:(_?BooleanQ):False] :=
Module[{strainEnergyDensity},
strainEnergyDensity=If[nodalField, MapToNodes[simData["strainEnergyDensity"], simData["modelData"]], simData["strainEnergyDensity"] ];
LocalUpdateRule[designData["designVector"], designData["penal"], designData["voidModulus"], strainEnergyDensity, utilizationMetric]
]

ToElemsIfNodal[nodalField_?BooleanQ, field_List, modelData_Association]:= If[nodalField, MapToElems[field, modelData], field]

(* overload SimulateFEModel from fem2d.wl to map design to nodes *)
SimulateFEModel[modelData_Association, designData_Association, nodalField:_?BooleanQ]:= Module[{mappedDesign},
mappedDesign = ToElemsIfNodal[nodalField, designData["designVector"], modelData];
simData = SimulateFEModel[modelData, UpdateKey[designData, "designVector", mappedDesign]]
]

(* assumption: the current design is implicitly defined in the designData data structure *)
DesignUpdate[modelData_Association, designData_Association, utilizationMetric_?NumericQ, nodalField:(_?BooleanQ):False] := Module[{simData, designVector},
simData = SimulateFEModel[modelData, designData, nodalField];
designVector = LocalUpdateRule[designData, simData, utilizationMetric, nodalField];
{designVector, simData}
]


(* TODO: if statement to set the numElems or numNodes as an Option, by default elems *)
InitialDesignVector[modelData_, initialValue_:1.]:= ConstantArray[initialValue, modelData["numElems"]] 

(*BUG: when OptimalityCriteria calls this function it always converges to all zeros!!!!*)
FixedPointSolver[residual_, initialGuess_, tol_:1*^-3]:= FixedPoint[residual, initialGuess, SameTest -> (Max[Abs[#1-#2]]<tol&)]



(* Volume Fraction Constraint *)

VolumeFractionResidual[designVector_List, volumeFraction_] := Mean[designVector] - volumeFraction

VolumeFractionResidual[designData_Association] := VolumeFractionResidual[designData["designVector"], designData["volumeFraction"]]

LagrangeMultiplierResidual[simData_Association, designData_Association, mu_, nodalField:(_?BooleanQ):False] := Module[{designVector},
designVector = LocalUpdateRule[designData, simData, mu, nodalField];
designVector = ToElemsIfNodal[nodalField, designVector, simData["modelData"]]; (* volume is computed using element volumes*)
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

DesignUpdate[modelData_Association, designData_Association, nodalField:(_?BooleanQ):False] := Module[{simData, mu, designVector},
simData = SimulateFEModel[modelData, designData, nodalField];
mu = BisectionRoot[LagrangeMultiplierResidual[simData, designData, #, nodalField]&][[1]];
designVector = LocalUpdateRule[designData, simData, mu, nodalField];
{designVector, simData, mu}
]




Options[OptimalityCriteria] = {"DensityFilter" -> Identity}
OptimalityCriteria[modelData_Association, designData_Association, OptionsPattern[]] := Module[{residual, elemField},
elemField = Length@designData["designVector"] == modelData["numElems"];
residual = OptionValue["DensityFilter"]@DesignUpdate[modelData, UpdateKey[designData, "designVector", #], !elemField][[1]]&;
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
