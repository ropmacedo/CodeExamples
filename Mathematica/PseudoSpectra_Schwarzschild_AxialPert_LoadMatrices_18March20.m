(* ::Package:: *)

(* ::Title:: *)
(*QNM Schwarzschild Potential*)


(* ::Section:: *)
(*Quasi-normal Modes and quasi-normal functions*)


(* ::Subsection:: *)
(*Numerical parameters*)


Clear["Global`*"]

SetDirectory[NotebookDirectory[]];
Prec=500;
z0=0;
z1=1;
\[CapitalDelta]z=z1-z0;
InputParFile=ToString[$CommandLine[[Length@$CommandLine]]];
InputFile=StringDelete[InputParFile,"Parameters/"];

Print[InputParFile];
Print[InputFile];
Get[InputParFile];



Nz=200;
spin=-2;
l=2;
Print[Nz];
Print[s];
NzHigh=2 Nz;
x[i_,Nz_]:=Cos[(\[Pi] i)/Nz];

X=N[Table[x[i,Nz],{i,0,Nz}],Prec];
XX=N[Table[x[i,NzHigh],{i,0,NzHigh}],Prec];

z=N[Table[z0+1/2 \[CapitalDelta]z (1+x[i,Nz]),{i,0,Nz}],Prec];
zz=N[Table[z0+1/2 \[CapitalDelta]z (1+x[i,NzHigh]),{i,0,NzHigh}],Prec];

\[Sigma]=z;(*(1+z)/2;*)
\[Sigma]\[Sigma]=zz;(*(1+zz)/2;*)


(* ::Subsubsection:: *)
(*Load Matrices*)


fn="OperatorMatrix/AxialParity/M_N_"<>ToString[Nz]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat"
M=N[Import[fn,"Table"]/10^(Prec+10),Prec];

fn="OperatorMatrix/AxialParity/MAdj0_N_"<>ToString[Nz]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat"
MAdj0=N[Import[fn,"Table"]/10^(Prec+10),Prec];


(* ::Subsubsection:: *)
(*Pseudo-Spectrum*)


dblId=IdentityMatrix[2*(Nz+1)];
A=M - s*dblId;

Adj0=MAdj0-Conjugate[s]*dblId;

Adj0A=Adj0 . A;



Print["Calculating Pseudo Spectra Energy Norm 0"];
PseudoSpectraAdj0A=Eigenvalues[Adj0A,-1];
PseudoSpectraAdj0AData={N@Re@s,N@Im@s,N[Abs@PseudoSpectraAdj0A,Prec]};
fn="Data/AxialParity/PseudoSpectraEnergyNorm0_spin_"<>ToString[spin]<>"_l_"<>ToString[l]<>"_N_"<>ToString[Nz]<>"_Prec_"<>ToString[Floor[Prec]]<>InputFile<>"_GnuPlot.dat";
Export[fn,PseudoSpectraAdj0AData];
Print["Done"];
