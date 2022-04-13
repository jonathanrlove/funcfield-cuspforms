//----------------------------------------------------------
//SECTIONS OF VECTOR BUNDLES
//----------------------------------------------------------
//
//

intrinsic MinimumEmptyToZero(seq::SeqEnum) -> RngIntElt
    {Returns the minimum of seq. If seq is empty, return 0.}

    if #seq eq 0 then
        return 0;
    else
        return Minimum(seq);
    end if;
end intrinsic;

intrinsic SectionSpace(adelmat::AdelMat : prec := 30) -> ., SeqEnum, SeqEnum, Map
    {Finds functions (f, h) which are holomorphic when multiplied by adelmat[v] for each v in S.
    
    To do: generalize to GL_n.}
    require Rank(adelmat) eq 2 : "Currently only functional for 2x2 adelic matrices.";
    C := Curve(adelmat);
    k := BaseRing(C);
    F := FunctionField(C);

    adelmat := Reduce(adelmat : scale := false);

    Df := DivisorGroup(C) ! 0;
    Dh := DivisorGroup(C) ! 0;
    for v in Support(adelmat) do
        Df +:= Valuation(Mat(adelmat, v)[1,1]) * Divisor(v);
        Dh +:= (Maximum(Valuation(Mat(adelmat, v)[1,1]) - Valuation(Mat(adelmat, v)[1,2]), 0) + Valuation(Mat(adelmat, v)[2,2])) * Divisor(v);
    end for;
    fspace, fmap := RiemannRochSpace(Df);
    hspace, hmap := RiemannRochSpace(Dh);
    fbasis := Basis(Df);
    hbasis := Basis(Dh);


    coeffspace := VectorSpace(k, #fbasis + #hbasis);

    for v in Support(adelmat) do
        Fv := Completion(v : prec:=prec);
        compmap := CompletionMap(v : prec:=prec);
        kv := ResidueClassField(v);
        conditionmat := Matrix(kv, #fbasis + #hbasis, 0, []);

        localmat := Mat(adelmat, v : prec:=prec);
        tv := localmat[1,2];
        alpha := Valuation(localmat[1,1]);
        beta := Valuation(tv);
        delta := Valuation(localmat[2,2]);

        minj := Minimum(MinimumEmptyToZero([Valuation(fi, v) + beta : fi in fbasis]), MinimumEmptyToZero([Valuation(hi, v) + delta: hi in hbasis]));

        for j in [minj..-1] do
            conditionvect := [];
            for fi in fbasis do
                fiexpans := compmap(fi);
                newterm := &+ChangeUniverse([Coefficient(fiexpans, l) * Coefficient(tv, j-l) : l in [Valuation(fiexpans) .. -beta-1]], kv);
                Append(~conditionvect, newterm);
            end for;
            for hi in hbasis do
                hiexpans := compmap(hi);
                Append(~conditionvect, Coefficient(hiexpans, j-delta));
            end for;
            conditionmat := HorizontalJoin(conditionmat, Transpose(Matrix(Vector(conditionvect))));
        end for;
        coeffspace := coeffspace meet RestrictField(Kernel(conditionmat), k);
    end for;

    toV := map<coeffspace -> VectorSpace(F, Rank(adelmat)) | 
               vect :-> [&+ChangeUniverse([vect[i] * fbasis[i] : i in [1..#fbasis]], F), //ChangeUniverse only necessary in case of summing over an empty list
                         &+ChangeUniverse([vect[#fbasis+i] * hbasis[i] : i in [1..#hbasis]], F)],
               Fvect :-> Eltseq(Fvect[1] @@ fmap) cat Eltseq(Fvect[2] @@ hmap)>;


    return coeffspace, toV;
end intrinsic;

intrinsic IsSection(sectionlist::SeqEnum, adelmat::AdelMat : prec := 30) -> BoolElt
    {Returns true iff every element of sectionlist is a section of adelmat.}
    C := Curve(adelmat);
    F := FunctionField(C);
    sectionmat := Matrix(sectionlist);
    require BaseRing(sectionmat) eq F : "funcpairlist must be defined over the function field of the base curve of adelmat";

    for w in Support(adelmat) join NonIntegralPlaces(sectionmat : checkdet := false) do
        localinfo := LocalizeMat(sectionmat, w : prec:=prec) * Mat(adelmat, w : prec:=prec);
        for f in ElementToSequence(localinfo) do
            if Valuation(f) lt 0 then
                return false;
            end if;
        end for;
    end for;
    return true;            
end intrinsic;

/*
intrinsic Section(vect::., fbasis::SeqEnum, hbasis::SeqEnum, adelmat::AdelMat : prec := 30) -> SeqEnum
    {vect should be a vector in the space output from SectionSpace}
    
    if #fbasis eq 0 then
        f := FunctionField(Curve(adelmat)) ! 0;
    else
        f := &+[vect[i] * fbasis[i] : i in [1..#fbasis]];
    end if;

    if #hbasis eq 0 then
        h := FunctionField(Curve(adelmat)) ! 0;
    else
        h := &+[vect[#fbasis+i] * hbasis[i] : i in [1..#hbasis]];
    end if;

    assert IsSection([[f,h]], adelmat : prec:=prec);

    return [f, h];
end intrinsic;
*/

intrinsic BasisOfSections(adelmat::AdelMat : prec := 30, oldsections := []) -> SeqEnum
    {Return a basis of global sections of adelmat. If oldsections is nonempty, return only as many sections as needed in order 
    to span the space of global sections over k together with oldsections.}
    F := FunctionField(Curve(adelmat));
    n := adelmat`Dim;

    C:= Curve(adelmat);
    space, toV := SectionSpace(adelmat : prec:=prec);

    oldpreims := [sect @@ toV : sect in oldsections];
    modifiedsections := ExtendBasis(oldpreims, space);
    newsections := [toV(vect) : vect in modifiedsections[#oldsections + 1..#modifiedsections]];

    //assert IsSection(newsections, adelmat : prec:=prec);
    return newsections;
end intrinsic;

//For isomorphism testing

intrinsic IndependentIndices(sectionlist::SeqEnum, admat::AdelMat) -> SeqEnum
    {A maximal list of indices i such that the set of all sectionlist[i] is independent.}
    F := FunctionField(Curve(admat));
    n := Rank(admat);
    require F eq BaseRing(Matrix(sectionlist)) : "sections must be defined over the function field of the base curve of admat";

    if #sectionlist eq 0 then
        return [], [];
    end if;

    V := VectorSpace(F, n);
    basis := [];
    indices := [];
    for i in [1..#sectionlist] do
        if IsIndependent(Append(basis, sectionlist[i])) then
            Append(~basis, sectionlist[i]);
            Append(~indices, i);
        end if;
    end for;

    return indices;
end intrinsic;

intrinsic ExtractCoeffs(mat::., i::RngIntElt) -> AlgMatElt
    {Given a matrix mat defined over a Laurent series ring with uniformizer pi, return the matrix of the coefficients of pi^i in each component.}

    Fv := BaseRing(mat);
    require Type(Fv) eq RngSerLaur : "mat must be defined over a Laurent Series Ring.";

    k := CoefficientRing(Fv);
    rows := NumberOfRows(mat);
    cols := NumberOfColumns(mat);
    return Matrix(k, rows, cols, [Coefficient(f, i) : f in ElementToSequence(mat)]);
end intrinsic;

intrinsic TwistAndCompare(admat1::AdelMat, admat2::AdelMat, twistamount::RngIntElt : prec := 30) -> BoolElt, AdelMat, AdelMat
    {Operates on admat1 and admat2 in place.}
    
    C := Curve(admat1);
    require C eq Curve(admat2) : "inputs must be defined over the same curve";
    v := Infinity(C);
    Fv := Completion(v : prec:=prec);

    SpecifyPlace(admat1, v, Mat(admat1, v) * Fv.1^twistamount);
    SpecifyPlace(admat2, v, Mat(admat2, v) * Fv.1^twistamount);
    dim1 := Dimension(SectionSpace(admat1 : prec:=prec));
    dim2 := Dimension(SectionSpace(admat2 : prec:=prec));
    return dim1, dim2;
end intrinsic;

intrinsic TwistForEnoughSections(admat1::AdelMat, admat2::AdelMat : prec := 30) -> BoolElt, AdelMat, AdelMat, SeqEnum, SeqEnum, SeqEnum
    {Test}

    require Curve(admat1) eq Curve(admat2) : "admat1 and admat2 must be defined over the same curve";
    require Rank(admat1) eq Rank(admat2) : "admat1 and admat2 must have the same rank";

    admat1 := Copy(admat1);
    admat2 := Copy(admat2);

    dim1 := Dimension(SectionSpace(admat1 : prec:=prec));
    dim2 := Dimension(SectionSpace(admat2 : prec:=prec));
    if dim1 ne dim2 then
        return false, admat1, admat2, [], [], [];
    end if;

    //Twist down until there are no sections
    while dim1 gt 0 do
        dim1, dim2 := TwistAndCompare(admat1, admat2, -1 : prec:=prec);
        if dim1 ne dim2 then
            return false, admat1, admat2, [], [], [];
        end if;
    end while;

    //Twist up until at least n+1 dimensions of sections over k and at least n dimensions of sections over F,
    //recording the number of dimensions obtained at each new step
    V := VectorSpace(admat1);
    n := Rank(admat1);
    basis1 := [];
    basis2 := [];
    bucketmax := [];
    while Dimension(sub<V | basis1>) lt n do
        newdim1, newdim2 := TwistAndCompare(admat1, admat2, 1 : prec:=prec);
        
        if newdim1 ne newdim2 then
            return false, admat1, admat2, [], [], [];
        end if;
        if newdim1 gt dim1 then
            newsections1 := BasisOfSections(admat1 : prec := prec, oldsections := basis1);
            newsections2 := BasisOfSections(admat2 : prec := prec, oldsections := basis2);
            bucketmax cat:= [#basis1 + #newsections1 : i in newsections1];
            basis1 cat:= newsections1;
            basis2 cat:= newsections2;
        end if;
        dim1 := newdim1;
        dim2 := newdim2;
    end while;

    return true, admat1, admat2, basis1, basis2, bucketmax;
end intrinsic;

intrinsic IsIsomorphic(admat1::AdelMat, admat2::AdelMat : prec:=30, allisos := false) -> BoolElt, SeqEnum
    {return true if admat1 and admat2 define isomorphic vector bundles, false if not. If true, also return the divisor D and
     matrix m in GL2(F) such that admat1 and m * admat2 * D are equivalent mod GL2(O) (ie up to "Reduce").}
    
    require Curve(admat1) cmpeq Curve(admat2) : "Inputs must be defined over the same curve";

    C := Curve(admat1);
    k := BaseRing(C);
    F := FunctionField(C);
    isolist := [];

    //If equivalent, then the difference of divisors is div(f)+2D
    Ddiff := DetDiv(admat1) - DetDiv(admat2);
    MultsofTwo, multbytwo := 2*Domain(Gtodiv(C));
    hashalf, originalhalf := HasPreimage(divtoG(C)(Ddiff), multbytwo);
    if not hashalf then
        return false, isolist;
    end if;

    newadmat1 := Reduce(admat1 : scale := false);
    redadmat1 := Reduce(admat1);

    for twotors in Kernel(multbytwo) do
        //modify admat2 so that the two adelic matrices have equal determinants
        half := originalhalf + twotors;
        newD := Gtodiv(C)(half);
        admat2withD := '*'(admat2, newD : Dsupp := ClassGroupRepSupport(C), prec := prec); //by supplying the support, we save on time that would be used calculating the support every runthrough
        tf, f := IsLinearlyEquivalent(DetDiv(admat1) - DetDiv(admat2), 2*newD);
        ymat := Matrix(F, 2, 2, [f, 0, 0, 1]);

        newadmat2 := Reduce('*'(ymat, admat2withD : prec := prec) : scale := false);
        assert DetDiv(newadmat1) eq DetDiv(newadmat2);

        samedims, twistadmat1, twistadmat2, basis1, basis2, bucketmax := TwistForEnoughSections(newadmat1, newadmat2 : prec := prec);

        if samedims then
            indepinds := IndependentIndices(basis1, twistadmat1);
            n := #indepinds;
            m := #basis1;

           
            /*
            if m eq n then bottomrows := Matrix(F, 0, n, []); 
            else bottomrows := Matrix([basis1[i] : i in [1..#basis1] | i notin indepinds]);
            end if;

            LeftMat := HorizontalJoin(bottomrows * toprows^-1, ScalarMatrix(F, m, -1));
            testvects := [];
            for i in [1..n+m] do
                for j in [1..n+m] do
                    testmat := ColumnSubmatrix(LeftMat, i, 1) * RowSubmatrix(Matrix(basis2), j, 1);
                    Append(~testvects, ElementToSequence(testmat));
                end for;
            end for;
            locmat := LocalizeMat(Matrix(testvects), Infinity(C) : prec:=prec);
            firstcoeff := Min([Valuation(f) : f in ElementToSequence(locmat)]);
            coeffsmat := ExtractCoeffs(locmat, firstcoeff);
            for i in [1..n+m+1] do
                coeffsmat := HorizontalJoin(coeffsmat, ExtractCoeffs(locmat, firstcoeff + i));
            end for;
            aoptions := Rows(NullspaceMatrix(coeffsmat));*/
            
            toprows := Matrix([basis1[i] : i in indepinds]);
            constraintmat := ZeroMatrix(k, n * m, 0);
            for i in [1..n] do
                for j in [bucketmax[indepinds[i]]+1..m] do
                    testmat := ZeroMatrix(k, n, m);
                    testmat[i,j] := 1;
                    constraintmat := HorizontalJoin(constraintmat, Matrix(k, n * m, 1, Eltseq(testmat)));
                end for;
            end for;
            aoptions := Rows(NullspaceMatrix(constraintmat));

            /*
            LeftMat := bottomrows * toprows^-1;
            ker := KernelMatrix(Matrix(basis2));
            locleftmat := LocalizeMat(LeftMat, Infinity(C) : prec:=prec);
            locker := LocalizeMat(ker, Infinity(C) : prec:=prec);
            firstcoeff := Min([Valuation(f) : f in ElementToSequence(locleftmat) cat ElementToSequence(locker)]);
            lastcoeff := Min([AbsolutePrecision(f) - 1 : f in ElementToSequence(locleftmat) cat ElementToSequence(locker)]);

            for d in [firstcoeff..lastcoeff] do
                if d ne 0 then
                    testvects := [];
                    left := ExtractCoeffs(locleftmat, d);
                    right := Transpose(KernelMatrix(Transpose(ExtractCoeffs(locker, d))));
                    for i in [1..n] do
                        for j in [1..bucketmax[n]] do
                            testmat := ColumnSubmatrix(left, i, 1) * RowSubmatrix(right, j, 1);
                            Append(~testvects, Eltseq(testmat));
                        end for;
                    end for;
                    constraintmat := HorizontalJoin(constraintmat, Matrix(testvects));
                end if;
            end for;
            aoptions1 := Rows(NullspaceMatrix(constraintmat));
            print aoptions1;
            */

            // Can we find a valid change of basis matrix?
            if #aoptions ge 1 then
                for coeffs in CartesianPower(k, #aoptions) do
                    samplesol := &+[coeffs[i] * aoptions[i] : i in [1..#aoptions]];
                    amat := Matrix(F, n, m, [samplesol[i] : i in [1..n*m]]);
                    changeofbasis := toprows^-1 * amat * Matrix(basis2);

                    if not IsSingular(changeofbasis) then
                        Fmat := changeofbasis * ymat;
                        tryadmat2 := Fmat * admat2withD;
                        
                        if AdMatsEqual(redadmat1, Reduce(tryadmat2)) then
                            Append(~isolist, <newD, Fmat>);
                            if not allisos then
                                return true, isolist;
                            end if;
                        end if;
                    end if; 
                end for;
            end if;        
        end if;
    end for;
    return #isolist gt 0, isolist;
end intrinsic;