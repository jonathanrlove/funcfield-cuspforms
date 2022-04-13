AttachSpec("Heckespec");

// Define a curve
Pol<x> := PolynomialRing(GF(3));
C := HyperellipticCurve(x^5+1,0);

// Define a list of places on C
placelist := RatPlaces(C); // includes all degree 1 places
F<x,y> := FunctionField(C);
Append(~placelist, Place([x^2-2*x+2,x-y]));

// loads sample Hecke graph information
load "SampleCurveData.m";
heckeinf := newhecke;

/* Computes Hecke graph information (type "HeckeSearch;" for details).
   WARNING: process takes several hours on a personal computer (use "tries"
   parameter to proceed in batches).
   Uncomment following lines if you're sure you want to use this. */

// heckeinf := NewHeckeInfo(C);
// HeckeSearch(heckeinf, placelist);

// Print information about simultaneous eigenspaces
ComputeEigenforms(heckeinf); // may take several minutes

load "InputOutput.m";
ToMathematica(heckeinf, "adjacencymatrix.nb" : Overwrite := true); //sends Hecke graphs to a Mathematica notebook
// OutputHeckeData(heckeinf, "SampleCurveData.m"); // Saves heckeinf to a file.
