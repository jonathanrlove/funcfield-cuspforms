//This package allows the user to work with curves (Crv) in a way that provides extra functionality.

declare attributes Crv: Gtodiv, divtoG, SpecialPlaces, Inf, Uniformizers, CompletionMaps, tests;

intrinsic Gtodiv(C::Crv) -> Map
    {A map from the (abstract) class group of C to the divisor group of C}
    if not assigned C`Gtodiv then
        G, Gtodiv, divtoG := ClassGroup(C);
        C`Gtodiv := Gtodiv;
        C`divtoG := divtoG;
    end if;
    return  C`Gtodiv;
end intrinsic

intrinsic divtoG(C::Crv) -> Map
    {A map from the divisor group of C to the (abstract) class group of C}
    if not assigned C`divtoG then
        G, Gtodiv, divtoG := ClassGroup(C);
        C`Gtodiv := Gtodiv;
        C`divtoG := divtoG;
    end if;
    return  C`divtoG;
end intrinsic

intrinsic ClassGroupRepSupport(C::Crv) -> SeqEnum
    {A list S of places of C such that every divisor on C is equivalent to one supported on S.}
    if not assigned C`SpecialPlaces then
        S := {};
        for g in Generators(Domain(Gtodiv(C))) do
            for P in Support(Gtodiv(C)(g)) do
                Include(~S, P);
            end for;
        end for;
        C`SpecialPlaces := S;
    end if;
    return C`SpecialPlaces;
end intrinsic;

intrinsic Infinity(C::Crv) -> PlcCrvElt
    {A pre-determined place of C.}
    if not assigned C`Inf then
        C`Inf := Random({v : v in ClassGroupRepSupport(C) | Degree(v) eq 1 });
    end if;
    return C`Inf;
end intrinsic;

intrinsic Uniformizer(v::PlcCrvElt) -> FldFunFracSchElt
    {Return a uniformizer u at v.}

    //Commented out:  If v is not in ClassGroupRepSupport(C), ensure div(u)-[v] is supported on ClassGroupRepSupport(C).
    
    C := Curve(v);

    if not assigned C`Uniformizers then
        C`Uniformizers := AssociativeArray(Places(C));
    end if;
    
    if not v in Keys(C`Uniformizers) then
        C`Uniformizers[v] := UniformizingParameter(v);
        /*
        if v in ClassGroupRepSupport(C) then
            C`Uniformizers[v] := UniformizingParameter(v);
        else
            vmovedtoS := Gtodiv(C)(divtoG(C)(Divisor(v)));
            tf, fnc := IsLinearlyEquivalent(Divisor(v), vmovedtoS);
            C`Uniformizers[v] := fnc;
        end if;
        */
    end if;
    return C`Uniformizers[v];
end intrinsic;

intrinsic CompletionMap(v::PlcCrvElt : prec := 30) -> Map
    {Return a map from the function field of C to its completion at v (with precision prec).
    The map is normalized so as to send Uniformizer(v) to the uniformizer of the completion.}

    C := Curve(v);

    if not assigned C`CompletionMaps then
        C`CompletionMaps := AssociativeArray(Places(C));
    end if;
    
    if not v in Keys(C`CompletionMaps) then
        F := FunctionField(C);
        Fv, phi := Completion(F, v : Precision:=prec);
        AssignNames(~Fv, ["t"]);

        piv := Uniformizer(v);
        intermsofpiv := Reversion(phi(piv));
        completionmap := map<F -> Fv | f :-> Composition(phi(f), intermsofpiv), locf :-> Composition(locf, phi(piv)) @@ phi>;
        C`CompletionMaps[v] := completionmap;
    end if;
    return C`CompletionMaps[v];
end intrinsic;

intrinsic Completion(v::PlcCrvElt : prec := 30) -> RngSerLaur
    {Return the completion of the function field of C at v (with precision prec).}
    return Codomain(CompletionMap(v : prec:=prec));
end intrinsic;

intrinsic UniformizedPlaces(C::Crv) -> SeqEnum
    {Returns a sequence of all places of C on which a uniformizer has been computed}
    pllist := [];
    if assigned C`Uniformizers then
        for v in Keys(C`Uniformizers) do
            Append(~pllist, v);
        end for;
    end if;
    return pllist;
end intrinsic;
/*
intrinsic LocalizeFunc(f::FldFunFracSchElt, v::PlcCrvElt : prec := 30) -> RngSerLaurElt
    {Returns the image of f (up to precision prec) in the completion of the function field of C at v. 
    The map is chosen such that Uniformizer(v) is sent to the uniformizer of the complete local ring.}
    require f in FunctionField(Curve(v)) : "f must be in the function field of C";
    piv := Uniformizer(v); 
    phi := CompletionMap(v : prec:=prec); 
    return Composition(phi(f), Reversion(phi(piv)));    
end intrinsic;

intrinsic PreimageFunc(f::RngSerLaurElt, v::PlcCrvElt : prec := 30) -> FldFunFracSchElt
    {Given f in the completion of a function field at v, 
    Returns an element of the function field reducing to f (up to precision prec).
    The map is chosen such that the uniformizer of the completion is sent to Uniformizer(v).}
    Fv := Completion(v : prec:=prec);
    require f in Fv : "f must be in the completion of the function field of C at v.";

    piv := Uniformizer(v); 
    phi := CompletionMap(v : prec:=prec); 
    return Composition(f, phi(piv)) @@ phi;
end intrinsic;
*/