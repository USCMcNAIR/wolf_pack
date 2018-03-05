(* ::Package:: *)

(* Modules for finite element analysis of rectangles of arbitrary aspect ratio and mesh size*)
(* by Luis Bahamonde Jacome *)
(* Copyright 2018  *)



BeginPackage["SquareFiniteElement`"]

FEM2D::usage = "FEM2D provides two interfaces for finite element analysis of rectangles of arbitrary aspect ratio and mesh size. 
                The high-level interface requires minimum understanding of the finite element analysis process and is documented in ?FEM2DHighLevel
                The low-level interface requires in-depth understanding of the finite element analysis process and is documented in ?FEM2DLowLevel"

FEM2DHighLevel::usage = " The high-level interface of the fem2d package provides data processing functions to model and simulate a parameterized rectangular domain.

The following data processing functions are available 

- modelData = ModelMBBBeam[inputData, designData]
- modelData = ModelCantileverBeam[inputData, designData]
- simData = SimulateFEModel[modelData, designData]

where the data structures are defined as associations with the following string keys:

- inputData: 'numElemsWidth', 'aspectRatio', 'youngModulus', 'poissonRatio'
- designData: 'voidModulus', 'penal', 'volumeFraction', 'designVector'
- modelData: intermediate variables necessary to solve for displacements
- simData: 'U', 'strainEnergy', 'strainEnergyDensity'

Note the modelData data structure contains as a key the inputData for traceability purposes. 
Note the simData data structure contains as a key the modelData for traceability purposes.
"

FEM2DLowLevel::usage = "The following low-level functions are available:
                         GetNumElems
                         GetNumDof
                         GenerateSquareElementMesh
                         ReshapeField
                         InterpolateYoungModulus
                         ElementLocalStiffnessMatrix
                         GlobalStiffnessVector
                         GlobalStiffnessMatrix
                         SolveDisplacements
                         StrainEnergyDensityField"
                         
GetNumElems::usage = "GetNumElems[numElemsWidth, aspectRatio] calculates the number of elements of the square mesh"

GetNumDof::usage = "GetNumDof[numElemsWidth, aspectRatio] calculates the total degrees of freedom in the mesh"

GenerateSquareElementMesh::usage = "GenerateSquareElementMesh[numElemsWidth, aspectRatio] creates a matrix of the global degrees of freedom that belong to each element. The element number corresponds to the row index of the matrix."

ReshapeField::usage = "ReshapeField[fieldVector, numElemsWidth] reshapes the field vector into a matrix with the same number of rows as numElemsWidth and the same number of columns as numElemsLength"

InterpolateYoungModulus::usage = "InterpolateYoungModulus[densityVector,  youngModulus, penal, voidModulus] interpolates between voidModulus and youngModulus based on the density vector and penalization"

ElementLocalStiffnessMatrix::usage = "ElementLocalStiffnessMatrix[nu] creates a bilinear quad element stiffness matrix, in local axes, normalized with respect to the Young modulus"
		
GlobalStiffnessVector::usage = "GlobalStiffnessVector[elementK, youngModulus, elementDofMatrix] creates a list of {posK, valK} for assembling the global stiffness matrix as a SparseArray.
                                Young modulus can be a scalar or a vector field of symbolically interpolated moduli.
                                If used for design, call it out of the optimization loop"	

GlobalStiffnessMatrix::usage = " GlobalStiffnessMatrix[numDof, posK, numK, densityVector] builds the global stiffness matrix as a numerical SparseArray. 
                                 If densityVector is given then numK must be a symbolic expression of \!\(\*SubscriptBox[\(\[Rho]\), \(i\)]\).
                                 If used for design, call it within the optimization loop, once the design variables of the current iteration are known" 
	
SolveDisplacements::usage = "SolveDisplacements[sparseK, fixedDof, sparseF] solves for the global displacement vector.
                             If sparseF is not given the output is the decomposed REDUCED stiffness matrix as a function to evaluate right-hand sides efficiently"

StrainEnergyDensityField::usage = "StrainEnergyDensityField[elementK, elementDofMatrix, U, youngModulus] provides a vector of the strain energy density of each element"

ModelRectangularStructure::usage = "ModelRectangularStructure[inputData] models a rectangular structure with no boundary condition or loading definitions"

ModelMBBBeam::usage = "ModelMBBBeam[inputData] models the classical MBB beam exploiting the symmetry of the problem"

ModelCantileverBeam::usage = "ModelCantileverBeam[inputData] models a cantilever beam with a load applied in the last degree of freedom, in a frowny-face bending sense"

SimulateFEModel::usage = "SimulateFEModel[modelData] solves an FE model data structure to obtain the displacements, and retrieve the strain energy"



Begin["`Private`"]


(* Preliminaries *)
SymbolicScalarField[scalar_, numElems_] := Subscript[scalar, #]& /@ Range[numElems]

DensityField[numElems_] :=(Clear[\[Rho]]; SymbolicScalarField[\[Rho], numElems])

DensityRules[densityField_List, designVector_List] := Dispatch@Thread[densityField -> designVector]

(* Element Mesh Generation *)

GetNumElemsLength[numElemsWidth_, aspectRatio_] := numElemsWidth * aspectRatio

GetNumElems[numElemsWidth_, aspectRatio_] := numElemsWidth * GetNumElemsLength[numElemsWidth, aspectRatio]

GetNumNodes[numElemsWidth_, aspectRatio_] := (numElemsWidth + 1) * (GetNumElemsLength[numElemsWidth, aspectRatio] + 1)

GetNumDof[numElemsWidth_, aspectRatio_] := 2 * GetNumNodes[numElemsWidth, aspectRatio]

GenerateSquareElementMesh[numElemsWidth_, aspectRatio_] := 
  Module[{nodesNumbering, edofVec},
  nodesNumbering = Partition[Range[GetNumNodes[numElemsWidth, aspectRatio]], numElemsWidth + 1]\[Transpose];
  edofVec = Flatten[(2 * nodesNumbering[[;; -2, ;;-2]]  + 1)\[Transpose]];
  Table[edofVec, 8]\[Transpose] + Table[Join[{0, 1}, 2 * numElemsWidth + {2, 3, 0, 1}, {-2, -1} ], GetNumElems[numElemsWidth, aspectRatio]]
]

ReshapeField[fieldVector_List, numElemsWidth_] := Transpose@Partition[fieldVector, numElemsWidth]


(* Material Properties *)
InterpolateYoungModulus[densityVector_List,  youngModulus_, penal_, voidModulus_] := voidModulus + Power[densityVector, penal] * (youngModulus - voidModulus)


(* FE Modeling *)
ElementLocalStiffnessMatrix[nu_Real] := Module[{A11, A12, B11, B12},
A11 = {{12, 3, -6, -3}, {3, 12, 3, 0}, {-6, 3, 12, -3}, {-3, 0, -3, 12}};
A12 = {{-6, -3, 0, 3}, {-3, -6, -3, -6}, {0, -3, -6, 3}, {3, -6, 3, -6}};
B11 = {{-4, 3, -2, 9}, {3, -4, -9, 4}, {-2, -9, -4, -3}, {9, 4, -3, -4}};
B12 = {{2, -3, 4, -9}, {-3, 2, 9, -2}, {4, 9, 2, 3}, {-9, -2, 3, 2}};
ArrayFlatten[(1/24)/(1 - nu^2) *({{A11, A12}, {Transpose[A12], A11}} + nu * {{B11, B12}, {Transpose[B12], B11}})]
]


(* FE Post processing *)
StrainEnergyDensity[elementK_List, elementU_List] := elementU.elementK.elementU 
 
StrainEnergyDensityField[elementK_List, elementDofMatrix_, U_, youngModulus_] := Module[{elementsUMatrix, evalSED},
elementsUMatrix = U[[#]]& /@ elementDofMatrix; 
evalSED = StrainEnergyDensity[elementK, #]&;
(evalSED /@ elementsUMatrix) * youngModulus
]


(* FE Analysis *)
GlobalStiffnessPositions[elementDofMatrix_List] := 
{ (*rowK*) Flatten[Transpose /@ Outer[Times, elementDofMatrix, ConstantArray[1, 8]], 1] //Flatten,
  (*colK*) Flatten /@ Outer[Times, elementDofMatrix, ConstantArray[1, 8]] //Flatten}

GlobalStiffnessVector[elementK_List, youngModulusField_List, elementDofMatrix_List] := Module[{rowK, colK, overlappingDof, valK},
{rowK, colK} = GlobalStiffnessPositions[elementDofMatrix];
overlappingDof = PositionIndex@Transpose[{rowK, colK}];
valK = Flatten@Transpose@Outer[Times, Flatten[elementK], youngModulusField];
{Keys[overlappingDof], Total[valK[[#]]]& /@ Values[overlappingDof]}]

GlobalStiffnessVector[elementK_List, youngModulus_Real, elementDofMatrix_List] := Module[{repeatedYoungModulus},
repeatedYoungModulus = ConstantArray[youngModulus, Length[elementDofMatrix]];
GlobalStiffnessVector[elementK, repeatedYoungModulus, elementDofMatrix]
]

GlobalStiffnessMatrix[numDof_, posK_List, numK_List] := Module[{sparseK},
sparseK = SparseArray[posK -> numK, {numDof, numDof}];
sparseK = 0.5 * (sparseK + Transpose[sparseK])
]

GlobalStiffnessMatrix[numDof_, posK_List, valK_List, densityRules_] := Module[{numK, sparseK},
numK = valK /. densityRules;
GlobalStiffnessMatrix[numDof, posK, numK]
]

GetFreeDof[fixedDof_, numDof_] := Complement[Range[numDof], fixedDof]

SolveDisplacements[sparseK_, fixedDof_] := Module[{freeDof},
freeDof = GetFreeDof[fixedDof, Dimensions[sparseK][[1]] (* numDof *)];
LinearSolve[sparseK[[freeDof, freeDof]], Method -> "Banded"]
]

SolveDisplacements[sparseK_, fixedDof_, sparseF_] := Module[{freeDof, numDof, U},
numDof = Dimensions[sparseK][[1]];
freeDof = GetFreeDof[fixedDof, numDof];
U = ConstantArray[0., numDof];
U[[freeDof]] = LinearSolve[sparseK[[freeDof, freeDof]], sparseF[[freeDof]], Method -> "Banded"];
U
]


(* Data processing*)

ReshapeField[simData_Association, fieldName_String] := ReshapeField[simData[fieldName], simData["modelData"]["inputData"]["numElemsWidth"]]

PrepareRectangularStructureModel[inputData_Association] := Module[{numElems, numDof, elementDofMatrix, elementK},
numElems = GetNumElems[inputData["numElemsWidth"], inputData["aspectRatio"]];
numDof = GetNumDof[inputData["numElemsWidth"], inputData["aspectRatio"]];
elementDofMatrix = GenerateSquareElementMesh[inputData["numElemsWidth"], inputData["aspectRatio"]];
elementK = ElementLocalStiffnessMatrix[inputData["poissonRatio"]];
AssociationThread[{"inputData", "numElems", "numDof", "elementDofMatrix", "elementK"} ->
                  { inputData,   numElems,   numDof,   elementDofMatrix,   elementK }]
]

ModelRectangularStructure[inputData_Association] := Module[{modelData, posK, numK, sparseK}, 
modelData = PrepareRectangularStructureModel[inputData];
{posK, numK} = GlobalStiffnessVector[modelData["elementK"], inputData["youngModulus"] , modelData["elementDofMatrix"]];
sparseK = GlobalStiffnessMatrix[modelData["numDof"], posK, numK];
AssociateTo[modelData, {"posK" -> posK, "numK" -> numK, "sparseK" -> sparseK}];
modelData
]

ModelRectangularStructure[inputData_Association, designData_Association] := Module[{modelData, densityField, youngModulusField, posK, valK, sparseK}, 
modelData = PrepareRectangularStructureModel[inputData];
densityField = DensityField[modelData["numElems"]];
youngModulusField = InterpolateYoungModulus[densityField, inputData["youngModulus"], designData["penal"], designData["voidModulus"]];
{posK, valK} = GlobalStiffnessVector[modelData["elementK"], youngModulusField , modelData["elementDofMatrix"]];
AssociateTo[modelData, {"designData" -> designData, "densityField" -> densityField, "youngModulusField" -> youngModulusField,
                        "posK" -> posK, "valK" -> valK, "sparseK" -> sparseK}];
modelData
]

ModelCantileverBeam[inputData_Association, designData_:0] := Module[{modelData, sparseF, fixedDof},
modelData = If[SameQ[designData, 0], ModelRectangularStructure[inputData], ModelRectangularStructure[inputData, designData]];
sparseF = SparseArray[{modelData["numDof"]} -> -1.0, {modelData["numDof"]}];
fixedDof = Range[1, 2 * (inputData["numElemsWidth"]+1), 1](* fixed end *);
AssociateTo[modelData, {"fixedDof" -> fixedDof, "sparseF" -> sparseF}];
modelData
]

ModelMBBBeam[inputData_Association, designData_:0] := Module[{modelData, sparseF, fixedDof},
modelData = If[SameQ[designData, 0], ModelRectangularStructure[inputData], ModelRectangularStructure[inputData, designData]];
sparseF = SparseArray[{2} -> -1.0, {modelData["numDof"]}];
fixedDof = Join[(* symmetry condition *) Range[1, 2*(inputData["numElemsWidth"] + 1), 2],(* edge pin *) {modelData["numDof"]}];
AssociateTo[modelData, {"fixedDof" -> fixedDof, "sparseF" -> sparseF}];
modelData
]

SimulateFEModel[modelData_Association]:= Module[{U, strainEnergyDensity, strainEnergy, simData},
U = SolveDisplacements[modelData["sparseK"], modelData["fixedDof"], modelData["sparseF"]];
strainEnergyDensity = StrainEnergyDensityField[modelData["elementK"], modelData["elementDofMatrix"], U, modelData["inputData", "youngModulus"]];
strainEnergy = Total[strainEnergyDensity];
simData = AssociationThread[{"modelData", "U", "strainEnergyDensity", "strainEnergy"} -> 
                            {modelData  ,  U ,  strainEnergyDensity ,  strainEnergy}];
simData
]

SimulateFEModel[modelData_Association, designData_Association] := Module[{densityRules, sparseK, U, interpolatedYoungModulus, strainEnergyDensity, strainEnergy, simData},
densityRules = DensityRules[modelData["densityField"], designData["designVector"]];
sparseK = GlobalStiffnessMatrix[modelData["numDof"], modelData["posK"], modelData["valK"], densityRules];
U = SolveDisplacements[sparseK, modelData["fixedDof"], modelData["sparseF"]];
interpolatedYoungModulus = modelData["youngModulusField"] /. densityRules;
strainEnergyDensity = StrainEnergyDensityField[modelData["elementK"], modelData["elementDofMatrix"], U, interpolatedYoungModulus];
strainEnergy = Total[strainEnergyDensity];
simData = AssociationThread[{"modelData", "designData", "densityRules",  "U", "strainEnergyDensity", "strainEnergy"} -> 
                            {modelData  ,  designData,   densityRules,    U ,  strainEnergyDensity ,  strainEnergy}];
simData
]


End[]
EndPackage[]

