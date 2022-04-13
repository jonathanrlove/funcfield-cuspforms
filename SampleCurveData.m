
//Data for curve:
C := HyperellipticCurve([Polynomial(GF(3), \[1, 0, 0, 0, 0, 1]), 
Polynomial(GF(3), \[])]);
FncField<X,Y> := FunctionField(C);
C`Inf := Place([X + 1, 2*Y]);

C`Uniformizers := AssociativeArray(Places(C));
C`Uniformizers[Place([X^3 + 2*X^2 + X + 1, Y + X + 1])] := 1/(X^2 + X)*Y + 1/X;
C`Uniformizers[Place([X^2 + X + 2, Y + 2*X])] := 1/(X^3 + 2*X^2 + 2*X + 2)*Y + 
2*X/(X^3 + 2*X^2 + 2*X + 2);
C`Uniformizers[Place([X^4 + 2*X^3 + X^2 + 2*X + 1, Y])] := 1/(X + 1)*Y;
C`Uniformizers[Place([X + 2, X^5 + 2])] := X + 2;
C`Uniformizers[Place([1/X, 2/X^3*Y])] := 2*X^2/(X^5 + 1)*Y;
C`Uniformizers[Place([X + 1, 2*Y])] := 2/(X^4 + 2*X^3 + X^2 + 2*X + 1)*Y;
C`Uniformizers[Place([X^2 + X + 2, Y + X])] := 1/(X^3 + 2*X^2 + 2*X + 2)*Y + 
X/(X^3 + 2*X^2 + 2*X + 2);
C`Uniformizers[Place([X, 2*Y + 1])] := 2/(X^4 + 2*X + 1)*Y + 1/(X^3 + 2*X^2 + X 
+ 1);
C`Uniformizers[Place([X^3 + 2*X^2 + X + 1, Y + 2*X + 2])] := 1/(X^2 + X)*Y + 
2/X;
C`Uniformizers[Place([X, 2*Y + 2])] := 2/(X^4 + 2*X + 2)*Y + (X + 2)/(X^4 + 2*X 
+ 2);

//Hecke Operator Data:
prec := 30;
newhecke := New(HeckeInfo);
newhecke`curve := C;
newhecke`spaces := 1000;
newhecke`index := 6;
newhecke`slants := [ 0, -1, -1, -1, -1, -2, 0, -2, 0, -2, -2, 2, -2, 0, 0, -2, 
4, -2, -2, 0 ];
newhecke`terminals := [];

R<t> := LaurentSeriesRing(BaseRing(C));

doubcosdata := [<[], []>,
<[Place([1/X, 2/X^3*Y]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([X, 2*Y + 1]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([X + 1, 2*Y]), Place([X, 2*Y + 2])], [Matrix(R, 2, 2, [t, 0, 0, t]), Matrix(R, 2, 2, [t, 0, 0, 1])]>,
<[Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t^2, 0, 0, t])]>,
<[Place([1/X, 2/X^3*Y]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t^2, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([1/X, 2/X^3*Y]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t^2, t, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([1/X, 2/X^3*Y]), Place([X, 2*Y + 1]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([1/X, 2/X^3*Y]), Place([X, 2*Y + 1]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t, 1, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([1/X, 2/X^3*Y]), Place([X + 1, 2*Y]), Place([X, 2*Y + 2])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t]), Matrix(R, 2, 2, [t, 0, 0, 1])]>,
<[Place([1/X, 2/X^3*Y]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t^2, 0, 0, t])]>,
<[Place([1/X, 2/X^3*Y]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t^2, t, 0, t])]>,
<[Place([X, 2*Y + 1]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t^2, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([X, 2*Y + 1]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t^2, t, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t])]>,
<[Place([X, 2*Y + 1]), Place([X + 1, 2*Y]), Place([X, 2*Y + 2])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t, 0, 0, t]), Matrix(R, 2, 2, [t, 1, 0, 1])]>,
<[Place([X, 2*Y + 1]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t^2, 0, 0, t])]>,
<[Place([X, 2*Y + 1]), Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t, 0, 0, 1]), Matrix(R, 2, 2, [t^2, t, 0, t])]>,
<[Place([X + 1, 2*Y]), Place([X, 2*Y + 2])], [Matrix(R, 2, 2, [t, 0, 0, t]), Matrix(R, 2, 2, [t^2, 0, 0, 1])]>,
<[Place([X + 1, 2*Y]), Place([X, 2*Y + 2])], [Matrix(R, 2, 2, [t^2, 0, 0, t]), Matrix(R, 2, 2, [t, 0, 0, 1])]>,
<[Place([X + 1, 2*Y])], [Matrix(R, 2, 2, [t^3, t^2, 0, t])]>];

newhecke`alldoubcos := [NewAdelMat(plcmat[1], plcmat[2], C : n := 2, prec := 
prec) : plcmat in doubcosdata];

newhecke`heckes := AssociativeArray(Places(C));
newhecke`heckes[Place([1/X, 2/X^3*Y])] := SparseMatrix(1000, 1000, [<1, 2, 4>, <2, 1, 1>, <2, 6, 1>, <2, 7, 2>, <3, 8, 1>, <3, 9, 3>, <4, 9, 3>, <4, 10, 1>, <5, 11, 1>, <5, 12, 3>]);
newhecke`heckes[Place([X, 2*Y + 1])] := SparseMatrix(1000, 1000, [<1, 3, 4>, <2, 8, 1>, <2, 9, 3>, <3, 1, 1>, <3, 13, 1>, <3, 14, 2>, <4, 6, 1>, <4, 15, 3>, <5, 16, 1>, <5, 17, 3>]);
newhecke`heckes[Place([X + 1, 2*Y])] := SparseMatrix(1000, 1000, [<1, 5, 4>, <2, 11, 1>, <2, 12, 3>, <3, 16, 1>, <3, 17, 3>, <4, 17, 3>, <4, 19, 1>, <5, 1, 1>, <5, 6, 1>, <5, 20, 2>]);
newhecke`heckes[Place([X, 2*Y + 2])] := SparseMatrix(1000, 1000, [<1, 4, 4>, <2, 9, 3>, <2, 10, 1>, <3, 6, 1>, <3, 15, 3>, <4, 1, 1>, <4, 14, 2>, <4, 18, 1>, <5, 17, 3>, <5, 19, 1>]);
