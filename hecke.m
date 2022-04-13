//Hecke operators

intrinsic SmallerSlant(adelmat::AdelMat : prec := 30) -> AdelMat
    {Return an equivalent adelic matrix with slant degree at most 2g.}

    C := Curve(adelmat);
    g := Genus(C);
    m := Floor(g - DegDet(adelmat)/2);
    F := FunctionField(C);

    if #Support(adelmat) eq 0 then
        return adelmat;
    end if;

    v := Infinity(C);
    Fv := Completion(v : prec := prec);
    SpecifyPlace(adelmat, v, Mat(adelmat, v) * Fv.1^m);
    assert DegDet(adelmat) in [2*g-1, 2*g];

    basis := BasisOfSections(adelmat : prec:=prec);
    assert #basis ge 1;
    
    nonzerof := false;
    for row in basis do
        if row[1] ne 0 then
            nonzerof := true;
            f := row[1];
            h := row[2];
            break;
        end if;
    end for;
    if nonzerof then
        adelmat := Matrix(F, 2, 2, [0, 1, f, h]) * adelmat;
    end if;

    assert DegSlant(Reduce(adelmat)) le 2*g;
    return Reduce(adelmat);
end intrinsic;

intrinsic HeckeNeighbors(adelmat::AdelMat, v::PlcCrvElt : prec := 30) -> SeqEnum
    {return the neighbors of adelmat according to the Hecke operator at v}
    require v in Places(Curve(adelmat)) : "v must be a place on the base curve of rep";
    C := Curve(adelmat);
    F := FunctionField(C);
    k := ResidueClassField(v);
    Fv := Completion(v : prec:=prec);

    leftcosets := [Matrix(Fv, 2,2,[Fv.1, i, 0, 1]) : i in k];
    Append(~leftcosets, Matrix(Fv, 2,2,[1,0,0,Fv.1]));

    neighbors := [];

    for leftcoset in leftcosets do
        neighbor := Copy(adelmat);
        SpecifyPlace(neighbor, v, Mat(neighbor, v) * leftcoset);
        Append(~neighbors, neighbor);
    end for;

    return neighbors;
    
end intrinsic;

intrinsic RatPlaces(C::Crv) -> SeqEnum
    {return a list of rational places of C.}

    return [Place(P) : P in RationalPoints(C)];
end intrinsic;

intrinsic HigherDegreePlace(C::Crv) -> PlcCrvElt
    {Returns a place of minimal degree greater than 1.}
    
    d := 2;
    while not HasPlace(C, d) do 
        d +:= 1; 
    end while;
    return RandomPlace(C, d);    
end intrinsic;

declare type HeckeInfo;

declare attributes HeckeInfo: curve, alldoubcos, allseenreps, slants, terminals, heckes, heckeplaces, spaces, eigenforms, eisensteins;

intrinsic NewHeckeInfo(C::Crv : spaces:=1000) -> HeckeInfo
    {Creates a new HeckeInfo object.}

    newinfo := New(HeckeInfo);
    newinfo`curve := C;
    newinfo`alldoubcos := [IdentityAdelMat(C)];
    newinfo`terminals := [];
    newinfo`slants := [0];
    newinfo`heckes := AssociativeArray(Places(C));
    newinfo`spaces := spaces;

    return newinfo;
end intrinsic;

intrinsic HeckeSearch(heckeinfo::HeckeInfo, places::SeqEnum : tries:=5000, prec:=30, terminalextension := 0, status := 2, checkback := false)
    {Takes in a HeckeInfo object over a curve C, and fills it with a list of double cosets on C
    and information about Hecke operators at each place in places. 
    
    More precisely, Heckesearch keeps a queue of all known double cosets. For each double coset 
    (at index i) in the queue, and each place v, the Hecke neighbors of adelmat at v are computed.
    They are then tested for isomorphism against all existing double cosets in the queue; if no
    isomorphisms are found, the neighbor is appended to the end of the queue. A matrix for each place
    is kept, where the entry at [i,j] is the number of Hecke neighbors from i to j.

    Optional parameters:
    - tries (default 5000): the number of pairs (i,v) (as above) to process. Note that if a pair
      has already been processed and stored in heckeinfo, then it will be skipped. This means that
      one can run HeckeSearch(heckeinfo, places : tries := 10); repeatedly, and each run will 
      process new information.

    - status (default 2) : an integer 1 (print i), 2 (print i and v), or 3 (print i, v, and all
      Hecke neighbors of i).
      
    Other parameters are not currently being used and should be left as default.}

    g := Genus(heckeinfo`curve);

    newplaces := [v : v in places | v notin Keys(heckeinfo`heckes)];
    for v in newplaces do
        heckeinfo`heckes[v] := SparseMatrix(heckeinfo`spaces, heckeinfo`spaces);
    end for;

    i := 1;
    Heckesprocessed := 0;
    while Heckesprocessed lt tries and i le #(heckeinfo`alldoubcos) do
        if status ge 1 then print "Double coset ",i,"/",#heckeinfo`alldoubcos," with slant =",heckeinfo`slants[i]; end if;
        for v in places do
            if NumberOfNonZeroEntries(RowSubmatrix(heckeinfo`heckes[v], i, 1)) eq 0 then // no Hecke neighbors computed yet
                Heckesprocessed +:= 1;
                if status ge 2 then print v; end if;
                for neighb in HeckeNeighbors(heckeinfo`alldoubcos[i], v : prec:=prec) do
                    new := true;
                    testj := [j : j in [1..i] | heckeinfo`heckes[v][j,i] ne 0]; // First check values of j that we know are sent to i by Hecke action.
                    testj cat:= [i..#heckeinfo`alldoubcos]; // Then check values of j beyond i. Hecke can send a double coset to itself!
                    testj cat:= [j : j in [1..i-1] | j in heckeinfo`terminals and heckeinfo`heckes[v][j,i] eq 0]; //We did not record all the hecke neighbors of terminals, so we should check those.
                    for j in testj do
                        tf, isolist := IsIsomorphic(neighb, heckeinfo`alldoubcos[j] : prec := prec);
                        if tf then
                            if status ge 3 then print j; end if;
                            heckeinfo`heckes[v][i,j] +:= 1;
                            new := false;
                            break;
                        end if;
                    end for;
                    if new then
                        slant := DegSlant(SmallerSlant(neighb : prec := prec));
                        if not i in heckeinfo`terminals or slant ge 2 - 2 * g - terminalextension then
                            Append(~heckeinfo`alldoubcos, neighb);
                            Append(~heckeinfo`slants, slant);
                            newindex := #heckeinfo`alldoubcos;
                            heckeinfo`heckes[v][i,newindex] +:= 1;
                            if status ge 3 then print "Found double coset", newindex, "with slant", slant; end if;
                            if slant lt 2 - 2 * g - terminalextension then
                                Append(~heckeinfo`terminals, newindex);
                            end if;
                        else
                            if status ge 3 then print "Found double coset with slant", slant,"adjacent to terminal"; end if;
                        end if;
                    end if;
                end for;
            end if;
        end for;
        i +:= 1;
    end while;
end intrinsic;

intrinsic HeckeMatrix(heckeinfo::HeckeInfo, pl::PlcCrvElt : sparse := false) -> AlgMatElt
    {Print the Hecke matrix at a given place}
    require pl in Keys(heckeinfo`heckes) : "heckeinfo must have been computed at place pl.";

    return Matrix(Submatrix(heckeinfo`heckes[pl], 1, 1, #heckeinfo`alldoubcos, #heckeinfo`alldoubcos));
end intrinsic;

intrinsic SimulEigens(space::ModTupRng, matlist::SeqEnum : status := 3) -> List
    {Computes simultaneous eigenvectors for a sequence of matrices acting on a space.
    If space can be decomposed into simultaneous eigenspaces, return true and a 
    sequence consisting of vectors and eigenvalues. Otherwise return false.}

    n := #matlist;
    
    if Dimension(space) eq 0 then
        return [* *];
    end if;

    if n eq 0 then
        return [* space *];
    end if;

    K := BaseRing(space);
    if status ge 1 then print "eigenspaces:", n, AbsoluteDegree(K); end if;
    mat := ChangeRing(matlist[1], K);
    proj := Matrix(Basis(space));
    projmat := proj * mat * Transpose(proj) * (proj * Transpose(proj))^-1;
    evalstocheck := GreatestCommonDivisor(CharacteristicPolynomial(mat), CharacteristicPolynomial(projmat));

    allsols := [* *];
    for factor in Factorization(evalstocheck) do
        fld := NumberField(factor[1]);
        eigenval := Roots(factor[1], fld)[1][1];

        newmat := ChangeRing(mat, fld);
        newspace := ChangeRing(space, fld);
        eigenspace := Eigenspace(newmat, eigenval) meet newspace;
        splitspace := SimulEigens(eigenspace, Remove(matlist, 1) : status := status);
        allsols cat:= splitspace;
    end for;
    
    return allsols;
end intrinsic;

intrinsic ComputeEigenforms(heckeinfo::HeckeInfo : status := 0, printeigens := true)
    {Compute eigenforms of heckeinfo, and assign them to heckeinfo`eigenforms.
    
    if status := 0 then no extra information is printed, if status is positive then
    updates are printed.

    If printeigens := true, print the eigenvalues at each place.}

    terminalmat := ZeroMatrix(IntegerRing(), #heckeinfo`alldoubcos, #heckeinfo`terminals);
    termindex := 1;
    for i in [1..#heckeinfo`terminals] do
        terminalmat[heckeinfo`terminals[i], i] := 1;
    end for;
    cuspformsupport := Kernel(terminalmat);
    
    mats := [Transpose(HeckeMatrix(heckeinfo, pl)) : pl in Keys(heckeinfo`heckes)]; // Magma computes left-eigenvectors rather than right

    heckeinfo`eigenforms := SimulEigens(cuspformsupport, mats : status := status);

    if printeigens then
        Pol<x> := PolynomialRing(Rationals());
        dims := [Dimension(eigenspace) : eigenspace in heckeinfo`eigenforms];
        if &*dims gt 1 then
            print "Simultaneous eigenspace dimensions too large: ", dims ,"Need Hecke operators at more places.";
        end if;
        for space in heckeinfo`eigenforms do
            for vect in Basis(space) do
                _, eigendata := IsSimultaneousEigenform(vect, heckeinfo);
                for v in Keys(eigendata) do
                    print v;
                    print AbsoluteMinimalPolynomial(eigendata[v]);
                end for;
            end for;
        end for;
    end if;
end intrinsic;

intrinsic IsSimultaneousEigenform(vect::ModTupFldElt, heckeinfo::HeckeInfo) -> BoolElt, SeqEnum, SeqEnum
    {Returns whether a form is an eigenform for Hecke operators at all places where heckeinfo has been computed.}

    fld := BaseRing(Parent(vect));
    data := AssociativeArray();
    for v in [v : v in Keys(heckeinfo`heckes)] do
        mat := Transpose(ChangeRing(HeckeMatrix(heckeinfo, v), fld));
        nonzeroi := [i : i in [1..#Eltseq(vect)] | vect[i] ne 0][1];
        eigenval := (vect * mat)[nonzeroi] / vect[nonzeroi];
        if IsZero(vect * mat - eigenval * vect) then
            data[v] := eigenval;
        else
            return false, _, _;
        end if;
    end for;
    return true, data;
end intrinsic;