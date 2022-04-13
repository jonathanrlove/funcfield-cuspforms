AttachSpec("Heckespec");

// Define a curve
Pol<x> := PolynomialRing(GF(3));
C := HyperellipticCurve(x^5+1,0);

// Define a list of places on C
placelist := RatPlaces(C); // includes all degree 1 places
F<x,y> := FunctionField(C);
Append(~placelist, Place([x^2-2*x+2,x-y]));

// Computes Hecke graph information (type ?HeckeSearch for details)
heckeinf := NewHeckeInfo(C);
HeckeSearch(heckeinf, placelist); 

// Print information about simultaneous eigenspaces
ComputeEigenforms(heckeinf); // also stores data in heckeinf`eigenforms

load "InputOutput.m";
ToMathematica(heckeinf, "adjacencymatrix.nb" : tag := "", Overwrite := false); //sends Hecke graphs to a Mathematica notebook


/* Other capabilities:

load "InputOutput.m";
OutputHeckeData(heckeinf, "SampleCurveData.m" : prec := 30) // Saves heckeinfo to a file.
load "SampleCurveData.m" // recovers heckeinfo, stored as variable "newhecke"
*/
