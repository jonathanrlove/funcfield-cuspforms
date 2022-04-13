174,0
S,LocalizeMat,Returns a matrix where FtoFv has been applied to every entry,0,2,0,0,0,0,0,0,0,232,,0,0,-34,,-34,-38,-38,-38,-38,-38
S,IsWeaklyEqual,Returns whether all corresponding components of mat1 and mat2 are weakly equal,0,2,0,0,0,0,0,0,0,177,,0,0,177,,36,-38,-38,-38,-38,-38
S,IsWeaklyIdentity,Returns whether mat is weakly equal to the identity matrix,0,1,0,0,0,0,0,0,0,177,,36,-38,-38,-38,-38,-38
S,NonIntegralPlaces,"Return the set of all places v where Fmat has non-integral components at v. If checkdet is true, also include all places v where Fmat has non-unit determinant (requires Fmat to be an invertible square matrix)",0,1,0,0,0,0,0,0,0,-34,,83,-38,-38,-38,-38,-38
S,EvaluateMat,Replace every instance of the uniformizer in mat with elt,0,2,0,0,0,0,0,0,0,285,,0,0,-34,,-34,-38,-38,-38,-38,-38
S,ReduceLocal,"Let Fv be a Laurent Series ring. The input mat must be an element of GL_2(Fv). Then there is r in GL_2(O_v) such that mat*r takes the form ((pi^a, p(pi)), (0, pi^d)) for integers a,d, and a Laurent polynomial p(X) in k[X] with no terms of degree a or higher. If scale=false, return mat*r, 1, and r. Additionally, there is a unique integer t such that mat*r*(pi^t) is integral and at least one component is a unit. If scale=true, return mat*r*(pi^t), pi^t, and r. To do: Generalize to GL_n(Fv)",0,1,0,0,0,0,0,0,0,177,,177,-38,-38,-38,-38,-38
T,AdelMat,-,0
A,AdelMat,3,Dict,BaseCurve,Dim
S,IdentityAdelMat,Return the identity adelic matrix on C,0,1,0,0,0,0,0,0,0,371,,AdelMat,-38,-38,-38,-38,-38
S,NewAdelMat,Creates a new adelic matrix over C with specified matrices at specified places,0,3,0,0,0,0,0,0,0,371,,0,0,82,,0,0,82,,AdelMat,-38,-38,-38,-38,-38
S,Curve,The base curve of adelmat,0,1,0,0,0,0,0,0,0,AdelMat,,371,-38,-38,-38,-38,-38
S,Rank,The rank (dimension) of adelmat,0,1,0,0,0,0,0,0,0,AdelMat,,148,-38,-38,-38,-38,-38
S,Dim,The dimension (rank) of adelmat,0,1,0,0,0,0,0,0,0,AdelMat,,148,-38,-38,-38,-38,-38
S,VectorSpace,"A vector space over the function field of the curve of adelmat, with specified dimension",0,1,0,0,0,0,0,0,0,AdelMat,,159,-38,-38,-38,-38,-38
S,Support,Return the set of places on which adelmat is nontrivial,0,1,0,0,0,0,0,0,0,AdelMat,,-50,-38,-38,-38,-38,-38
S,Mat,Returns the matrix at place v of the adelic matrix adelmat,0,2,0,0,0,0,0,0,0,232,,0,0,AdelMat,,177,-38,-38,-38,-38,-38
S,AdMatsEqual,Returns whether admat1 and admat2 are weakly equal at every place,0,2,0,0,0,0,0,0,0,AdelMat,,0,0,AdelMat,,36,-38,-38,-38,-38,-38
S,eq,Returns whether admat1 and admat2 are weakly equal at every place,0,2,0,0,0,0,0,0,0,AdelMat,,0,0,AdelMat,,36,-38,-38,-38,-38,-38
S,Print,"Print admat at level L. If L is ""Minimal,"" the dimensions, base curve, and support of admat are printed. Otherwise, the matrix at each place in the support of admat is also printed",0,2,0,0,1,0,0,0,0,298,,0,0,AdelMat,,-38,-38,-38,-38,-38,-38
S,AdelMatData,Prints the structures needed to define adelmat using the AdelMat intrinsic. Only works if all places in the support of adelmat are rational,0,1,0,0,0,0,0,0,0,AdelMat,,82,82,371,148,-38,-38
S,Copy,return a distinct adelic matrix that is identical to adelmat (so that changing the return value does not change adelmat),0,1,0,0,0,0,0,0,0,AdelMat,,AdelMat,-38,-38,-38,-38,-38
S,SpecifyPlace,Changes adelmat to set the matrix at v equal to mat,0,3,0,0,1,0,0,0,0,177,,0,0,232,,0,0,AdelMat,,-38,-38,-38,-38,-38,-38
S,*,Returns the product of the two adelic matrices,0,2,0,0,0,0,0,0,0,AdelMat,,0,0,AdelMat,,AdelMat,-38,-38,-38,-38,-38
S,*,Returns the product of a matrix over the function field by an adelic matrix,0,2,0,0,0,0,0,0,0,AdelMat,,0,0,-34,,AdelMat,-38,-38,-38,-38,-38
S,*,Returns the product of an adelic matrix by a divisor,0,2,0,0,0,0,0,0,0,61,,0,0,AdelMat,,AdelMat,-38,-38,-38,-38,-38
S,Reduce,"At every place v, right-multiply adelmat by an element of GL2(O_v) such that the result is upper triangular, with powers of the uniformizer on the diagonals. If scale = true, also multiply by an element of F_v to ensure that all components are integral, and at least one is a unit. Return the result",0,1,0,0,0,0,0,0,0,AdelMat,,AdelMat,-38,-38,-38,-38,-38
S,SlantDiv,"Reduces adelmat, then returns the divisor defined at each place to be the valuation of component (1,1) minus the valuation of component (2,2)",0,1,0,0,0,0,0,0,0,AdelMat,,61,-38,-38,-38,-38,-38
S,DegSlant,The degree of SlantDiv(adelmat),0,1,0,0,0,0,0,0,0,AdelMat,,148,-38,-38,-38,-38,-38
S,DetDiv,The divisor defined at each place to be the valuation of the determinant,0,1,0,0,0,0,0,0,0,AdelMat,,61,-38,-38,-38,-38,-38
S,DegDet,The degree of DetDiv(adelmat),0,1,0,0,0,0,0,0,0,AdelMat,,148,-38,-38,-38,-38,-38
S,EisensteinOld,Compute Eisenstein series,0,1,0,0,0,0,0,0,0,AdelMat,,-1,-38,-38,-38,-38,-38
S,Eisenstein,Compute Eisenstein series,0,2,0,0,0,0,0,0,0,172,,0,0,AdelMat,,-1,-38,-38,-38,-38,-38
S,Chi,compute character,0,2,0,0,0,0,0,0,0,172,,0,0,AdelMat,,-1,-38,-38,-38,-38,-38
S,NaiveEisenstein,Compute Eisenstein series,0,2,0,0,0,0,0,0,0,172,,0,0,AdelMat,,-1,-38,-38,-38,-38,-38
