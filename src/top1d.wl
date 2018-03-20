(* ::Package:: *)

(* MODULES for Truss Topology OPtimization *)
(* by Luis Bahamonde *)
(* Copyright 2016  *)
(* McNair Aerospace Research Center *)
(* University of South Carolina *)
(* e-mail: luisb@email.sc.edu *) 
 
BeginPackage["TrussOptimTop`"]

TrussOptimTop::usage = 
"  The following modules are available in the Truss Topology Opotimization package:
    ** showGroundStructure
** OptimalTruss
"

showGroundStructure::usage = 
" showGroundStructure[nx, ny, distx, disty, DX, DY]  shows the 2-D ground structure made from bars and nodes:
nx     -- number of nodes along x
ny     -- number of nodes along y
distx  -- neghbor connections along x
disty  -- neghbor connections along y
DX     -- node stepsize along x
DY     -- node stepzise along y
"

OptimalTruss::usage = 
"
OptimalTruss[xmax, ymax, nx, ny, supports, loads, \[IndentingNewLine]                    distx_:  1, disty_:0, kappa_: 1,  tol_: Sqrt[$MachineEpsilon]]

xmax -- rectangular design domain length 
ymax -- rectangular design domain width
nx -- number of nodes along x
ny -- number of nodes along y
supports -- displacement boundary conditions: {spprt1, spprt2, ...} = {{{node_id},{constraint_switch}}, ...} = {{{id_x, id_y},{ux, uy}}, ...}
loads -- applied loads:  {F1, F2, F3, ...} = {{{application_point}, {force_vector}}, ...} =  {  {{x1, y1},{Nx1, Ny1}}, {{x2, y2},{Nx2, Ny2}}, {{x3, y3},{Nx3, Ny3}}, ...  }
"


Needs["Developer`"] ;

(* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =*)
(* - - - TRUSS NODE AND BAR MESH - - - - - - - - - - - - - - - - - *)
(* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =*)


GenInc[dist_] := ToPackedArray[ Flatten[ Reap[
Sow[  { {1, 0},{0, 1},{1, 1},{-1, 1} }  ];
Do[   If[ GCD[i,j] == 1,
   Sow[{{i,j}, {-i,j}, {j,i}, {-j,i}}]]
      , {i, 2, dist}, {j, 1, i-1}  ]] [[2]], 2]]
(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -*)

Patterns[nx_, ny_, distx_, disty_] := Module[
{idx, idy, i0, i1, ne, INC},
INC = GenInc[Max[distx, disty]];
ToPackedArray[Reap[Do[
{idx, idy} = INC[[i]];
{i0, i1} = If[idx>= 0, {0, nx - idx},{-idx, nx}];
If[i0 <= i1 && idy <= ny && 
  Abs[idx] <= distx && idy <= disty, 
     Sow[{idx, idy, i0, i1, (i1 - i0 +1)  (ny - idy + 1)}]]
        , {i, Length[INC]}]][[2,1]], Integer]]
(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -*)

ElXY [nx_, ny_, patt_, DX_, DY_] := Module[
{idx, idy, i0, i1, ne, X, Y},
X = ToPackedArray[Table[i DX, {i,0, nx}], Real];
Y = ToPackedArray[Table[j DY, {j, 0, ny}], Real];
ToPackedArray[Flatten[Reap[
Do[{idx, idy, i0, i1, ne} = patt[[p]];
Sow[     Table[ {{X[[i]], Y[[j]]}, {X[[i +idx]], Y[[j + idy]] }}
 , {j, 1, ny - idy + 1},{i, i0 +1, i1 + 1} ]     ];
    , {p, Length[patt]}]][[2]], 3], Real]] 

(* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -*)
showGroundStructure[nx_, ny_, distx_, disty_, DX_, DY_] :=
		Graphics[ Line[     ElXY[nx, ny, Patterns[nx, ny, distx, disty ], DX, DY]    ]]




(* Lengths and direction cosines *)
LCS[patt_, DX_, DY_] := Module[{L,c,s},
ToPackedArray[Table[
 c = patt[[p, 1]]  DX; s = patt[[p,2]]  DY;
L = Sqrt[c*c + s*s]; c/= L; s /= L;
{L, c, s}
 , {p, Length[patt]}], Real]]


VectorL[patt_, lcs_] := Module[{ne, L},
ToPackedArray[ Flatten [ Table[
ne = patt[[p, 5]]; L = lcs[[p,1]]; Table[L,{ne}]
, {p, Length[patt]}]], Real]]



(* Boundary Conditions *)
BCList[supports_, nx_] := Module[{k, i, j, ux, uy}, 
ToPackedArray[ Partition[Union[ Flatten[ Reap[
Do[ {{i,j},{ux, uy}} = supports[[s]];
k = 2 ((nx + 1) j + i) + 1;
If[ux == 1, Sow[k]]; If[uy == 1, Sow[k+1]];
, {s, Length[supports]}]] [[2,1]] ]], 1], Integer]]


(* direction cosines B assembly*)
MatrixBT[nx_, ny_, patt_, lcs_] := Module[
{idx, idy, i0, i1, L, c, s, rules, k, dk, ie = 0},
rules = Reap[Do[{idx, idy, i0, i1, k} = patt[[p]];
{L, c, s} = lcs[[p]];
dk = 2 ((nx + 1) idy + idx);
Do[ie++; k = 2 ((nx + 1) j + i)  + 1;
Sow[ {ie,k} -> -c]; Sow[{ie, k+1} -> -s];
k += dk;
Sow[ {ie,k} -> c]; Sow[{ie, k+1} -> s];
 , {j,0, ny - idy},{i, i0, i1}];
  , {p, Length[patt]}];] [[2,1]];
Transpose[SparseArray[rules]]]

MatrixBT[nx_, ny_, patt_, lcs_, BC_] :=
Delete[MatrixBT[nx, ny, patt, lcs], BC]

(* Nodal Forces P *)
VectorP[loads_, nx_, ny_] := Module[{P, i, j, n, p},
P = Table[{0,0},{(nx + 1) (ny + 1)}];
Do[{{i, j},p} = loads [[f]];
    n = (nx + 1) j + i + 1; P[[n]] = p;
    , {f, Length[loads]}];
ToPackedArray[Flatten[P], Real]]

VectorP[loads_, nx_, ny_, BC_] := ToPackedArray[
Delete[VectorP[loads, nx, ny], BC], Real]


(* Topology Optimization *)
OptimalTruss[xmax_, ymax_, nx_, ny_, supports_, loads_, 
                    distx_:  1, disty_:0, kappa_: 1,  tol_: Sqrt[$MachineEpsilon]]:= 
Module[{ndx, ndy, dx, dy, pat, lcs, L, LL, BC, BB, PP,
              nn, ne, dof, S, A, e, g, a, t, vol, P, B},

ndx = Min[nx, Max[1, distx]];
ndy = Min[ny, Max[1, If[disty < 1, distx, disty]]];
dx   =   xmax/nx; dy = ymax / ny ;
pat = Patterns[nx, ny, ndx, ndy];
nn = (nx + 1)(ny + 1);
ne = Total[   pat[[All,5]]  ];
BC = BCList[supports, nx];
dof = 2 nn - Length[BC] ;
Print["Mesh ", nx, "x", ny, ":", ndx, "x", ndy,
", Nodes ", nn, ", Elements ", ne, ", DOF ", dof];
lcs = LCS[pat, dx, dy];
L = VectorL[pat, lcs] ;
P = VectorP[loads, nx, ny, BC];
B = MatrixBT[nx, ny, pat, lcs, BC] ;

PP = Transpose[ {   P, Table[0, { Length[P] }]          }];
LL = Join[L, kappa L];
BB = Join [B, -B, 2];
Print["Matrix H ", Length[P], "x", Length[LL], " in ",
  ByteCount[BB] / 2.^20, "MB (", 16  dof  ne / 2.^30, "GB full)" ];
t= Timing[S = LinearProgramming[LL, BB, PP, 
   Method -> "InteriorPoint", Tolerance->tol];] [[1]];
vol = S.LL;
 Print[ "Objective S.L = ", vol, " CPU time = ", t, "s"];
S = Partition[S, ne];
S = S[[1]] - S[[2]];
A = Table[  If[S[[i]] < 0, -kappa S[[i]], S[[i]] ] , {i, ne}  ];
A/= Max[A];
A = Chop[A, 100 tol];

e = ElXY[nx, ny, pat, dx, dy];
Graphics[ Reap[ Do[a = A[[i]]; If[a > tol,
    Sow[{Thickness[.015 Sqrt[a]], Hue[.7 (1 - a)],
Line[  e[[i]]  ]}]], {i,ne}]]  [[2, 1]]  ]       ]



EndPackage[]

