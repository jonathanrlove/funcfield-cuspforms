// Procedures for handling input/output

procedure PrintList(string, data : breaks := 100)
    //prints data according to formatting given by string.

    printf "[";
    done := 0;
    for elt in data do
        if done ne 0 and done mod breaks eq 0 then
            printf ",\n    ";
        elif done ne 0 then
            printf ", ";
        end if;
        printf string, elt;
        done +:= 1;
    end for;
    printf "]";
end procedure;

function TwoGensList(v)
    //Returns a pair of generators.

    f, g := TwoGenerators(v);
    return "[" cat Sprint(f) cat ", " cat Sprint(g) cat "]";
end function;

procedure OutputHeckeData(heckeinfo, file : prec := 30, startcomment := "", Overwrite := true)
    //Prints the info necessary to duplicate heckeinfo

    Write(file, startcomment : Overwrite := Overwrite);
    SetOutputFile(file);
    
    print "//Data for curve:";
    printf "C := %m;\n", heckeinfo`curve;

    F := FunctionField(heckeinfo`curve);
    AssignNames(~F, ["X","Y"]);
    printf "FncField<X,Y> := FunctionField(C);\n";
    printf "C`Inf := Place(%o);\n\n", TwoGensList(heckeinfo`curve`Inf);

    printf "C`Uniformizers := AssociativeArray(Places(C));\n";
    for v in Keys(heckeinfo`curve`Uniformizers) do
        printf "C`Uniformizers[Place(%o)] := %o;\n", TwoGensList(v), heckeinfo`curve`Uniformizers[v];
    end for;

    print "\n//Hecke Operator Data:";
    printf "prec := %o;\n", prec;
    printf "newhecke := New(HeckeInfo);\n";
    printf "newhecke`curve := C;\n";
    printf "newhecke`spaces := %o;\n", heckeinfo`spaces;
    printf "newhecke`slants := %o;\n", heckeinfo`slants;
    printf "newhecke`terminals := %o;\n\n", heckeinfo`terminals;
    printf "R<t> := LaurentSeriesRing(BaseRing(C));\n\n";
    

    first := true;
    printf "doubcosdata := [<";
    for admat in heckeinfo`alldoubcos do
        
        if not first then
            printf "]>,\n<";
        end if;
        first := false;

        
        plclist, matlist := AdelMatData(admat);
        //SetOutputFile(file);
        PrintList("Place(%o)", [TwoGensList(v) : v in plclist]);
        printf ", [";

        innerfirst := true;
        for mat in matlist do
            if not innerfirst then
                printf ", ";
            end if;
            innerfirst := false;
            printf "Matrix(R, 2, 2, ";
            PrintList("%m", Eltseq(mat));
            printf ")";
        end for;
    end for;
    printf "]>];\n\n";

    printf "newhecke`alldoubcos := [NewAdelMat(plcmat[1], plcmat[2], C : n := 2, prec := prec) : plcmat in doubcosdata];\n\n";

    //UnsetOutputFile(); //In case an error occurs
    printf "newhecke`heckes := AssociativeArray(Places(C));\n";
    for v in Keys(heckeinfo`heckes) do
        printf "newhecke`heckes[Place(%o)] := SparseMatrix(%o, %o, ", TwoGensList(v), 1000, 1000; //heckeinfo`spaces, heckeinfo`spaces;
        PrintList("%o", Eltseq(heckeinfo`heckes[v]));
        printf ");\n";
    end for;

    /*
    if assigned heckeinfo`eigenforms then
        printf "\n\n\nnewhecke`eigenforms := [* *];\n";
        for eigenspace in heckeinfo`eigenforms do
            printf "lst := [";
            first := true;
            for vect in Basis(eigenspace) do
                if first then first := false; 
                else printf ", "; end if;
                printf "%m", Eltseq(vect);
            end for;
            printf "];\nAppend(~newhecke`eigenforms, sub<VectorSpace(Universe(lst[1]), #(lst[1])) | lst>);\n\n";
        end for;
    end if;
    */

    UnsetOutputFile();
end procedure;

procedure ToMathematica(heckeinfo, file : tag := "", Overwrite := true)
    //Outputs the Hecke matrices to be readable by Mathematica

    Write(file, "" : Overwrite := Overwrite);
    SetOutputFile(file);

    printf "slants%o = {", tag;
    first := true;
    for elt in heckeinfo`slants do 
        if not first then printf ","; end if; first := false;
        printf "%o", elt;
    end for;
    printf "};\n\n";

    counter := 1;
    for v in Keys(heckeinfo`heckes) do
        printf "Mat%opl%o = SparseArray[{", tag, counter;
        first := true;
        for triad in Eltseq(heckeinfo`heckes[v]) do
            if not first then printf ", "; end if; first := false;
            printf "{%o,%o}->%o", triad[1], triad[2], triad[3];
        end for;
        printf "}, {%o,%o}];\n\n", #heckeinfo`alldoubcos, #heckeinfo`alldoubcos;
        printf "AdjacencyGraph[Mat%opl%o, VertexLabels -> Table[i -> slants%o[[i]], {i, 1, %o}]]\n\n", tag, counter, tag, #heckeinfo`alldoubcos;
        counter +:= 1;
    end for;

    printf "\n\n\n";

    UnsetOutputFile();
end procedure;

procedure WriteEigenData(heckeinf, eigendatafile)
    //Write Eigendata

    SetOutputFile(eigendatafile);
    printf "//New Curve\n\n";
    printf "C := %m;\n", heckeinf`curve;

    F := FunctionField(heckeinf`curve);
    AssignNames(~F, ["X","Y"]);
    printf "FncField<X,Y> := FunctionField(C);\n";
    printf "Append(~Curves,C);\n";

    places := [v : v in Keys(heckeinf`heckes)];
    printf "PlaceLst := ";
    PrintList("Place(%o)", [TwoGensList(v) : v in places]);

    printf ";\neigendata := [* *];\n\n";
    first := true;
    for eigenspace in heckeinf`eigenforms do
        _, eigenlist, vlist := IsSimultaneousEigenform(Basis(eigenspace)[1], heckeinf);
        neweigenlist := [eigenlist[Index(vlist, v)] : v in places];
        printf "Append(~eigendata, <%o, %o, %m>);\n", Dimension(eigenspace), AbsoluteDegree(Universe(neweigenlist)), neweigenlist;
    end for;
    printf "Append(~CurveData, <PlaceLst, eigendata>);\n\n\n", heckeinf`curve;
    UnsetOutputFile();
end procedure;

procedure ProcessCurves(Clist, eigendatafile : start := 1, prec := 30)
    //Compute Hecke operators and eigenforms

    if start eq 1 then 
        Write(eigendatafile, "" : Overwrite := true);
        SetOutputFile(eigendatafile);
        printf "Curves := [];\n";
        printf "CurveData := [* *];\n\n";
        UnsetOutputFile();
    end if;

    heckesearchtime := []; // [ 8272.547, 7157.813, 8480.156, 10989.312 ]
    eigenformstime := []; // [ 42.078, 33.672, 4.359, 2.766 ]
    eigenvaldata := [];
    for i in [start..#Clist] do
        C := Clist[i];
        clgroupsize := #TorsionSubgroup(ClassGroup(C));
        newprec := Maximum(prec, 2 * clgroupsize + 10);
        heckeinf := NewHeckeInfo(C);

        tim := Cputime();
        HeckeSearch(heckeinf, RatPlaces(C) : tries := 1000, status := 1, prec := newprec);
        v := HigherDegreePlace(C);
        print "New place of degree", Degree(v);
        upto := heckeinf`index;
        heckeinf`index := 1;
        HeckeSearch(heckeinf, [v] : tries := upto - 1, status := 1, prec := newprec);
        HeckeSearch(heckeinf, RatPlaces(C) cat [v] : tries := 1000, status := 1, prec := newprec);
        tim1 := Cputime(tim);
        Append(~heckesearchtime, tim1);

        tim := Cputime();
        ComputeEigenforms(heckeinf : status := 1);
        tim2 := Cputime(tim);
        Append(~eigenformstime, tim2);
        
        comment := "//Hecke search: " cat Sprint(tim1) cat " seconds\n";
        comment cat:= "//Eigenforms: " cat Sprint(tim2) cat " seconds\n\n";

        OutputHeckeData(heckeinf, "Curve" cat Sprint(i) cat "Data.m" : startcomment := comment);
        WriteEigenData(heckeinf, eigendatafile);
        
        print heckesearchtime;
        print eigenformstime;
    end for;
end procedure;
