//This package sets up the basic functionality of adelic matrices.

//The next several functions give functionality for matrices over the function field of a curve, matrices over a 
//completion, and methods for going between the two.

intrinsic LocalizeMat(mat::Mtrx, v::PlcCrvElt : prec := 30) -> Mtrx
    {Returns a matrix where FtoFv has been applied to every entry.}
    require BaseRing(mat) cmpeq FunctionField(Curve(v)) : "components of mat must be in the function field of C";

    Fv := Completion(v : prec:=prec);
    C := Curve(v);
    lst := [CompletionMap(v : prec:=prec)(f) : f in ElementToSequence(mat)];
    return Matrix(Fv, Nrows(mat), Ncols(mat), lst);
end intrinsic;

intrinsic IsWeaklyEqual(mat1::AlgMatElt, mat2::AlgMatElt) -> BoolElt
    {Returns whether all corresponding components of mat1 and mat2 are weakly equal.}
    require Nrows(mat1) eq Nrows(mat2) : "Dimensions are not equal";
    require Ncols(mat1) eq Ncols(mat2) : "Dimensions are not equal";
    for i in [1..Nrows(mat1)] do
        for j in [1..Ncols(mat1)] do
            if not IsWeaklyEqual(mat1[i,j], mat2[i,j]) then
                return false;
            end if;
        end for;
    end for;
    return true;
end intrinsic;

intrinsic IsWeaklyIdentity(mat::AlgMatElt) -> BoolElt
    {Returns whether mat is weakly equal to the identity matrix.}
    if Ncols(mat) ne Nrows(mat) then
        return false;
    end if;

    return IsWeaklyEqual(mat, ScalarMatrix(BaseRing(mat), Nrows(mat), 1));
end intrinsic;

intrinsic NonIntegralPlaces(Fmat::Mtrx : checkdet := true) -> SetEnum
    {Return the set of all places v where Fmat has non-integral components at v. If checkdet is true, 
    also include all places v where Fmat has non-unit determinant (requires Fmat to be an invertible square matrix).}
    places := [];
    for f in ElementToSequence(Fmat) do
        if f ne 0 then
            places cat:= Poles(f);
        end if;
    end for;
    if checkdet then
        require Nrows(Fmat) eq Ncols(Fmat) : "Fmat must be a square matrix";
        det := Determinant(Fmat);
        require det ne 0 : "Fmat must be invertible";
        places cat:=Poles(det);
        places cat:=Zeros(det);
    end if;
    return Set(places);
end intrinsic;

intrinsic EvaluateMat(mat::Mtrx, elt::RngSerLaurElt) -> Mtrx
    {Replace every instance of the uniformizer in mat with elt.}

    require Type(BaseRing(mat)) eq RngSerLaur : "mat must be defined over a Laurent Series Ring.";

    newmatelts := [Evaluate(f, elt) : f in ElementToSequence(mat)];
    return Matrix(Parent(elt), Nrows(mat), Ncols(mat), newmatelts);
end intrinsic;

intrinsic ReduceLocal(mat::AlgMatElt : scale := true) -> AlgMatElt
    {Let Fv be a Laurent Series ring. The input mat must be an element of GL_2(Fv). Then there is r in GL_2(O_v)
    such that mat*r takes the form ((pi^a, p(pi)), (0, pi^d)) for integers a,d, and a Laurent polynomial p(X) in k[X]
    with no terms of degree a or higher. If scale=false, return mat*r, 1, and r.

    Additionally, there is a unique integer t such that mat*r*(pi^t) is integral and at least one component is a unit. 
    If scale=true, return mat*r*(pi^t), pi^t, and r.

    To do: Generalize to GL_n(Fv).}

    require Nrows(mat) eq 2 : "mat must be a 2x2 matrix";
    require Ncols(mat) eq 2 : "mat must be a 2x2 matrix";
    require not IsWeaklyEqual(Determinant(mat), 0) : "mat must be invertible";
    Fv := BaseRing(mat);
    require Type(Fv) eq RngSerLaur : "mat must be defined over a Laurent series ring";

    newmat := mat;
    Kmat := Matrix(Fv, 2, 2, [1, 0, 0, 1]);
    //Ensure that newmat[2,2] has valuation no larger than newmat[2,1]
    if Valuation(newmat[2,1]) lt Valuation(newmat[2,2]) then
        newmat *:= Matrix(Fv, 2, 2, [0, 1, -1, 0]);
        Kmat *:= Matrix(Fv, 2, 2, [0, 1, -1, 0]);
    end if;

    //Make the matrix triangular
    Kmat *:= Matrix(Fv, 2, 2, [1, 0, -newmat[2,1]/newmat[2,2], 1]);
    newmat *:= Matrix(Fv, 2, 2, [1, 0, -newmat[2,1]/newmat[2,2], 1]);

    //Make the diagonals pure powers of Fv.1
    oneoverap := Fv.1 ^ (Valuation(newmat[1,1])) / newmat[1,1];
    oneoverdp := Fv.1 ^ (Valuation(newmat[2,2])) / newmat[2,2];
    Kmat *:= Matrix(Fv, 2, 2, [oneoverap, 0, 0, oneoverdp]);
    newmat *:= Matrix(Fv, 2, 2, [oneoverap, 0, 0, oneoverdp]);

    //Reduce newmat[1,2] mod newmat[1,1]
    bp := (newmat[1,2] - Truncate(newmat[1,2] + BigO(newmat[1,1]))) / Fv.1^Valuation(newmat[1,1]);
    Kmat *:= Matrix(Fv, 2, 2, [1, -bp, 0, 1]);
    newmat *:= Matrix(Fv, 2, 2, [1, -bp, 0, 1]);

    t := 0;
    if scale then
        //Scale
        t := -Valuation(newmat[2,2]); //To make bottom right component equal to 1
        //t := -Min([Valuation(newmat[1,1]), Valuation(newmat[1,2]), Valuation(newmat[2,2])]); //To make it integral
        newmat *:= Fv.1^t;
    end if;

    //newmat is a product of mat, t, and r:
    assert IsWeaklyEqual(mat * Fv.1^t * Kmat, newmat);
    //r=Kmat is in GL_2(O_v):
    assert Valuation(Kmat[1,1]) ge 0 and Valuation(Kmat[1,2]) ge 0 and Valuation(Kmat[2,1]) ge 0 and Valuation(Kmat[2,2]) ge 0;
    assert Valuation(Kmat[1,1]*Kmat[2,2] - Kmat[1,2]*Kmat[2,1]) eq 0;
    //newmat has the desired form (up to known precision):
    assert IsWeaklyEqual(newmat, Matrix(Fv, 2, 2, [Fv.1^Valuation(newmat[1,1]), Truncate(newmat[1,2] + BigO(newmat[1,1])), 0, Fv.1^Valuation(newmat[2,2])]));

    return newmat, Fv.1^t, Kmat;
end intrinsic;

//----------------------------------------------------------
//ADELIC MATRICES
//----------------------------------------------------------
//
//The type AdelMat allows the user to work with finitely supported elements of GL_n(A_F), where F is the function field
//of a projective curve and A_F is the adeles of F. Such elements can be represented by a finite set of places v,
//together with a matrix in GL_n(F_v) for each v in this set. The element is assumed to equal the identity matrix at 
//all other places of the curve.
//
//Note: the code does not require elements to be invertible.

declare type AdelMat;

declare attributes AdelMat: Dict, BaseCurve, Dim;

intrinsic IdentityAdelMat(C::Crv : n := 2, prec := 30) -> AdelMat
    {Return the identity adelic matrix on C.}

    adelmat := New(AdelMat);
    adelmat`BaseCurve := C;
    adelmat`Dim := n;
    adelmat`Dict := AssociativeArray(Places(C));

    return adelmat;
end intrinsic;

intrinsic NewAdelMat(plclist::SeqEnum, matlist::SeqEnum, C::Crv : n := 2, prec := 30) -> AdelMat
    {Creates a new adelic matrix over C with specified matrices at specified places.}

    require #plclist eq #matlist : "Must be as many places as matrices.";
    adelmat := New(AdelMat);
    adelmat`BaseCurve := C;
    adelmat`Dim := n;
    adelmat`Dict := AssociativeArray(Places(C));

    if #plclist eq 0 then
        return adelmat;
    end if;

    require Universe(plclist) eq Places(C) : "plclist must be a list of places on C.";
    require Type(Universe(matlist)) eq AlgMat : "matlist must be a list of matrices";

    for ind in [1..#plclist] do
        v := plclist[ind];
        locmat := matlist[ind];
        require Type(BaseRing(locmat)) eq RngSerLaur : "Matrices must be defined over a Laurent Series Ring";
        require Nrows(locmat) eq n : "Dimension mismatch (input n must equal dimensions of matrices)";
        require Ncols(locmat) eq n : "Dimension mismatch (input n must equal dimensions of matrices)";

        Fv := Completion(v : prec:=prec);
        adelmat`Dict[v] := EvaluateMat(locmat, Fv.1); 
    end for;

    return adelmat;
end intrinsic;

intrinsic Curve(adelmat::AdelMat) -> Crv
    {The base curve of adelmat}
    return adelmat`BaseCurve;
end intrinsic;

intrinsic Rank(adelmat::AdelMat) -> RngIntElt
    {The rank (dimension) of adelmat}
    return adelmat`Dim;
end intrinsic;

intrinsic Dim(adelmat::AdelMat) -> RngIntElt
    {The dimension (rank) of adelmat}
    return adelmat`Dim;
end intrinsic;

intrinsic VectorSpace(adelmat::AdelMat) -> ModTupFld
    {A vector space over the function field of the curve of adelmat, with specified dimension}
    return VectorSpace(FunctionField(adelmat`BaseCurve), adelmat`Dim);
end intrinsic;

intrinsic Support(adelmat::AdelMat) -> Set
    {Return the set of places on which adelmat is nontrivial.}
    return Keys(adelmat`Dict);
end intrinsic;

intrinsic Mat(adelmat::AdelMat, v::PlcCrvElt : prec := 30) -> AlgMatElt
    {Returns the matrix at place v of the adelic matrix adelmat}
    require v in Places(Curve(adelmat)) : "v must be a place of the base curve of adelmat";
    if v in Support(adelmat) then
        return adelmat`Dict[v];
    else
        Fv := Completion(v : prec := prec);
        return ScalarMatrix(Fv, adelmat`Dim, 1);
    end if;
end intrinsic;

intrinsic AdMatsEqual(admat1::AdelMat, admat2::AdelMat) -> BoolElt
    {Returns whether admat1 and admat2 are weakly equal at every place.}
    require admat1`BaseCurve cmpeq admat2`BaseCurve : "Inputs must be defined over the same curve";
    require admat1`Dim cmpeq admat2`Dim : "Dimensions must be equal"; 
    
    if Support(admat1) ne Support(admat2) then
        return false;
    end if;

    for v in Support(admat1) do
        if not IsWeaklyEqual(Mat(admat1, v), Mat(admat2, v)) then
            return false;
        end if;
    end for;

    return true;
end intrinsic;

intrinsic 'eq'(admat1::AdelMat, admat2::AdelMat) -> BoolElt
    {Returns whether admat1 and admat2 are weakly equal at every place.}

    return AdMatsEqual(admat1, admat2);
end intrinsic;

intrinsic Print(admat::AdelMat, L::MonStgElt)
    {Print admat at level L. If L is "Minimal," the dimensions, base curve, and support of admat are printed.
    Otherwise, the matrix at each place in the support of admat is also printed.}
    
    if #Support(admat) eq 0 then
        printf "%o x %o identity adelic matrix over curve %o", admat`Dim, admat`Dim, admat`BaseCurve;
    else
        if L eq "Minimal" then
            printf "%o x %o adelic matrix over %o supported at %o", admat`Dim, admat`Dim, admat`BaseCurve, Support(admat);
        else
            printf "%o x %o adelic matrix over %o with decomposition:\n", admat`Dim, admat`Dim, admat`BaseCurve;
            for v in Support(admat) do
                print v, admat`Dict[v];
            end for;
        end if;
    end if;
end intrinsic;

intrinsic AdelMatData(adelmat::AdelMat) -> SeqEnum, SeqEnum, Crv, RngIntElt
    {Prints the structures needed to define adelmat using the AdelMat intrinsic.
    Only works if all places in the support of adelmat are rational.}

    plclist := [];
    matlist := [* *];
    C := adelmat`BaseCurve;
    n := adelmat`Dim;

    for v in Support(adelmat) do
        Append(~plclist, v);
        R<t> := LaurentSeriesRing(ResidueClassField(v));
        Append(~matlist, EvaluateMat(adelmat`Dict[v], t));
    end for;

    return plclist, matlist, C, n;
end intrinsic;

intrinsic Copy(adelmat::AdelMat) -> AdelMat
    {return a distinct adelic matrix that is identical to adelmat (so that changing the return value does not change adelmat)}

    newadmat := New(AdelMat);
    newadmat`BaseCurve := adelmat`BaseCurve;
    newadmat`Dim := adelmat`Dim;
    newadmat`Dict := adelmat`Dict;

    return newadmat;
end intrinsic;

intrinsic SpecifyPlace(adelmat::AdelMat, v::PlcCrvElt, mat::AlgMatElt)
    {Changes adelmat to set the matrix at v equal to mat.}
    require v in Places(Curve(adelmat)) : "v must be a place on the base curve of adelmat";
    require BaseRing(mat) cmpeq Completion(v) : "mat must be defined over the completion at v of the function field of the base ring of adelmat";

    adelmat`Dict[v] := mat;
end intrinsic;

intrinsic '*'(admat1::AdelMat, admat2::AdelMat : prec := 30) -> AdelMat
    {Returns the product of the two adelic matrices.}
    require Curve(admat1) cmpeq Curve(admat2) : "Factors must be defined over the same curve";
    require admat1`Dim cmpeq admat2`Dim : "Dimensions not compatible for matrix multiplication";

    newadmat := New(AdelMat);
    newadmat`BaseCurve := admat1`BaseCurve;
    newadmat`Dim := admat1`Dim;
    newadmat`Dict := AssociativeArray(Places(newadmat`BaseCurve));

    for v in Support(admat1) join Support(admat2) do
        newadmat`Dict[v] := Mat(admat1, v : prec:=prec) * Mat(admat2, v : prec:=prec);
    end for;

    return newadmat;
end intrinsic;

intrinsic '*'(Fmat::Mtrx, admat2::AdelMat : prec := 30) -> AdelMat
    {Returns the product of a matrix over the function field by an adelic matrix.}
    //require BaseRing(Fmat) cmpeq FunctionField(Curve(admat2)) : "Factors must be defined over the same curve";
    //require Ncols(Fmat) cmpeq admat2`Dim : "Dimensions not compatible";
    //require Ncols(Fmat) cmpeq Nrows(Fmat) : "Fmat must be a square matrix";
    //require not IsSingular(Fmat) : "Fmat must be invertible";

    newadmat := New(AdelMat);
    newadmat`BaseCurve := Curve(admat2);
    newadmat`Dim := admat2`Dim;
    newadmat`Dict := AssociativeArray(Places(Curve(admat2)));

    placestocheck := Support(admat2) join NonIntegralPlaces(Fmat : checkdet := true);
    
    for v in placestocheck do
        newadmat`Dict[v] := LocalizeMat(Fmat, v : prec := prec) * Mat(admat2, v : prec:=prec);
    end for;

    return newadmat;
end intrinsic;

intrinsic '*'(admat1::AdelMat, D::DivCrvElt : prec := 30, Dsupp := {}) -> AdelMat
    {Returns the product of an adelic matrix by a divisor.}

    newadmat := New(AdelMat);
    C := Curve(admat1);
    newadmat`BaseCurve := C;
    newadmat`Dim := admat1`Dim;
    newadmat`Dict := AssociativeArray(Places(C));

    if #Dsupp eq 0 then Dsupp := Set(Support(D)); end if;

    for v in Support(admat1) join Dsupp do
        Fv := Completion(v : prec := prec);
        newadmat`Dict[v] := Mat(admat1, v : prec:=prec) * Fv.1^(Valuation(D, v));
    end for;

    return newadmat;
end intrinsic;

intrinsic Reduce(adelmat::AdelMat : scale := true) -> AdelMat
    {At every place v, right-multiply adelmat by an element of GL2(O_v) such that the result is
    upper triangular, with powers of the uniformizer on the diagonals. If scale = true, also multiply
    by an element of F_v to ensure that all components are integral, and at least one is a unit.
    Return the result.}

    newadmat := Copy(adelmat);
    for v in Support(adelmat) do
        newadmat`Dict[v] := ReduceLocal(adelmat`Dict[v] : scale := scale);
        if IsWeaklyIdentity(newadmat`Dict[v]) then
            Remove(~newadmat`Dict, v);
        end if;
    end for;

    return newadmat;
end intrinsic;

intrinsic SlantDiv(adelmat::AdelMat) -> DivCrvElt
    {Reduces adelmat, then returns the divisor defined at each place to be the valuation of component (1,1) minus
    the valuation of component (2,2).}
    require adelmat`Dim eq 2 : "adelmat must be a 2 x 2 adelic matrix.";

    adelmat := Reduce(adelmat);
    D := DivisorGroup(Curve(adelmat)) ! 0;
    for v in Support(adelmat) do
        D +:= (Valuation(adelmat`Dict[v][1,1]) - Valuation(adelmat`Dict[v][2,2])) * Divisor(v);
    end for;
    return D;
end intrinsic;

intrinsic DegSlant(adelmat::AdelMat) -> RngIntElt
    {The degree of SlantDiv(adelmat).}
    return Degree(SlantDiv(adelmat));
end intrinsic;

intrinsic DetDiv(adelmat::AdelMat) -> DivCrvElt
    {The divisor defined at each place to be the valuation of the determinant.}

    D := DivisorGroup(Curve(adelmat)) ! 0;
    for v in Support(adelmat) do
        D +:= Valuation(Determinant(Mat(adelmat, v))) * Divisor(v);
    end for;
    return D;
end intrinsic;

intrinsic DegDet(adelmat::AdelMat) -> RngIntElt
    {The degree of DetDiv(adelmat).}
    return Degree(DetDiv(adelmat));
end intrinsic;

intrinsic EisensteinOld(adelmat::AdelMat) -> .
    {Compute Eisenstein series}

    C := Curve(adelmat);
    todiff := FunctionField(C).1; // to compute a differential

    Ga, isom := AdditiveGroup(GF(3)); // defining a character
    psi := CharacterTable(Ga)[2];
    Cycfld := CoefficientField(psi);
    Cycpol<var> := PolynomialRing(Cycfld);
    Cycfldext<sqrt> := ext<Cycfld | var^2 - 3>;

    twos := 4;
    
    admat := Reduce(adelmat);
    K := Divisor(Differential(todiff));
    RR, modtofunc := RiemannRochSpace(K + DetDiv(admat)); // functions such that Whittaker is nonzero

    constant1 := 1;
    constant2 := Evaluate(ZetaFunction(C), 1/3^(twos-1)) / Evaluate(ZetaFunction(C), 1/3^(twos));
    for v in Support(admat) do
        constant1 *:= sqrt^(- Degree(v) * twos * Valuation(Mat(admat, v)[1,1]));
        constant2 *:= sqrt^(- Degree(v) * (2 - twos) * Valuation(Mat(admat, v)[1,1]));
    end for;

    total := constant1 + constant2;

    for i in RR do
        xi := modtofunc(i);
        if xi ne 0 then
            term := 1;
            for v in Support(admat) join Set(Support(Divisor(xi))) join Set(Support(K)) do
                locxi := CompletionMap(v)(xi);
                omega := Derivative(CompletionMap(v)(todiff));
                x := Mat(admat, v)[1,2];
                y := Mat(admat, v)[1,1];
                q := #ResidueClassField(v);

                mv := Valuation(locxi*y*omega);
                res := AbsoluteTrace(Coefficient(locxi*x*omega, -1));

                if mv lt 0 then
                    term := 0;
                else
                    term *:= psi(res @@ isom) * (sqrt^(-Degree(v) * (2 - twos) * Valuation(y))) * (1 - q^((mv + 1) * (1 - twos))) / (1 - q^(1 - twos));
                end if;
            end for;
            total +:= term / Evaluate(ZetaFunction(C), 1/3^(twos));
        end if;
    end for;
    return total;
/*
    total := 0;
    for i in RR do
        xi := modtofunc(i);
        if xi ne 0 then
            term := 1;
            for v in Support(admat) join Set(Support(Divisor(xi))) join Set(Support(K)) do
                locxi := CompletionMap(v)(xi);
                omega := Derivative(CompletionMap(v)(todiff));
                x := Mat(admat, v)[1,2];
                y := Mat(admat, v)[1,1];

                mv := Valuation(locxi*y*omega);
                res := AbsoluteTrace(Coefficient(locxi*x*omega, -1));

                term *:= psi(res @@ isom) * (sqrt^(-Degree(v)*mv)) * Max(0,mv+1);
            end for;
            total +:= term;
        end if;
    end for;
    return total;
    */
end intrinsic;

intrinsic Eisenstein(adelmat::AdelMat, s::FldComElt) -> .
    {Compute Eisenstein series}

    C := Curve(adelmat);
    todiff := FunctionField(C).1; // to compute a differential

    Ga, isom := AdditiveGroup(GF(3)); // defining a character
    psi := CharacterTable(Ga)[2];
    
    admat := Reduce(adelmat);
    K := Divisor(Differential(todiff));
    //RR, modtofunc := RiemannRochSpace(DetDiv(admat)); // functions such that Whittaker is nonzero
    RR, modtofunc := RiemannRochSpace(K + DetDiv(admat)); // functions such that Whittaker is nonzero

    constant1 := 3^(Genus(C) - 1);
    constant2 :=  Evaluate(ZetaFunction(C), 1/3^(2*s-1)) / Evaluate(ZetaFunction(C), 1/3^(2*s));
    for v in Support(admat) do
        constant1 *:= (#ResidueClassField(v)) ^ ((- s) * Valuation(Mat(admat, v)[1,1]));
        constant2 *:= (#ResidueClassField(v)) ^ ((s - 1) * Valuation(Mat(admat, v)[1,1]));
    end for;

    total := constant1 + constant2;

    for i in RR do
        xi := modtofunc(i);
        if xi ne 0 then
            term := 1;
            for v in Support(admat) join Set(Support(Divisor(xi))) join Set(Support(K)) do
                locxi := CompletionMap(v)(xi);
                omega := Derivative(CompletionMap(v)(todiff));
                x := Mat(admat, v)[1,2];
                y := Mat(admat, v)[1,1];
                q := #ResidueClassField(v);

                mv := Valuation(locxi*y*omega);
                res := AbsoluteTrace(Coefficient(locxi*x*omega, -1));

                if mv lt 0 then
                    term := 0;
                else
                    term *:= Conjugates(psi(res @@ isom))[1] * q^((s - 1) * Valuation(y)) * (1 - q^((mv + 1) * (1 - 2*s))) / (1 - q^(1 - 2*s));
                end if;
            end for;
            //print Abs(term / Evaluate(ZetaFunction(C), 1/3^(2*s))), ",";
            total +:= term / Evaluate(ZetaFunction(C), 1/3^(2*s));
        end if;
    end for;
    return total;
/*
    total := 0;
    for i in RR do
        xi := modtofunc(i);
        if xi ne 0 then
            term := 1;
            for v in Support(admat) join Set(Support(Divisor(xi))) join Set(Support(K)) do
                locxi := CompletionMap(v)(xi);
                omega := Derivative(CompletionMap(v)(todiff));
                x := Mat(admat, v)[1,2];
                y := Mat(admat, v)[1,1];

                mv := Valuation(locxi*y*omega);
                res := AbsoluteTrace(Coefficient(locxi*x*omega, -1));

                term *:= psi(res @@ isom) * (sqrt^(-Degree(v)*mv)) * Max(0,mv+1);
            end for;
            total +:= term;
        end if;
    end for;
    return total;
    */
end intrinsic;

intrinsic Chi(adelmat::AdelMat, s::FldComElt) -> .
    {compute character}

    admat := Reduce(adelmat);
    total := 1;
    for v in Support(admat) do
        total *:= (#ResidueClassField(v)) ^ (- s * Valuation(Mat(admat, v)[1,1]));
        total *:= (#ResidueClassField(v)) ^ (s * Valuation(Mat(admat, v)[2,2]));
    end for;
    return total;
end intrinsic;

intrinsic NaiveEisenstein(adelmat::AdelMat, s::FldComElt) -> .
    {Compute Eisenstein series}

    C := Curve(adelmat);
    F<X> := FunctionField(C);
    Pol<x> := PolynomialRing(GF(3));

    total := Chi(adelmat, s);
    for r in GF(3) do
        Fmat := Matrix(F,2,2,[0,-1,1,r]);
        total +:= Chi(Fmat * adelmat, s);
        print total;
    end for;
    for c in CartesianPower(GF(3), 3) do
        num := c[1]*x+c[2];
        den := c[3]*x+1;
        if num ne 0 and GCD(num, den) eq 1 then
            r := Evaluate(num, X) / Evaluate(den, X);
            print r;
            Fmat := Matrix(F,2,2,[0,-1,1,r]);
            total +:= Chi(Fmat * adelmat, s);
            print total;
            r := 1/r;
            Fmat := Matrix(F,2,2,[0,-1,1,r]);
            total +:= Chi(Fmat * adelmat, s);
            print total;
        end if;
    end for;
    return total;
end intrinsic;