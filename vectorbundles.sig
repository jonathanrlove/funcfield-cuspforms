174,0
S,MinimumEmptyToZero,"Returns the minimum of seq. If seq is empty, return 0",0,1,0,0,0,0,0,0,0,82,,148,-38,-38,-38,-38,-38
S,SectionSpace,"Finds functions (f, h) which are holomorphic when multiplied by adelmat[v] for each v in S. To do: generalize to GL_n",0,1,0,0,0,0,0,0,0,AdelMat,,-1,82,82,175,-38,-38
S,IsSection,Returns true iff every element of sectionlist is a section of adelmat,0,2,0,0,0,0,0,0,0,AdelMat,,0,0,82,,36,-38,-38,-38,-38,-38
S,BasisOfSections,"Return a basis of global sections of adelmat. If oldsections is nonempty, return only as many sections as needed in order to span the space of global sections over k together with oldsections",0,1,0,0,0,0,0,0,0,AdelMat,,82,-38,-38,-38,-38,-38
S,IndependentIndices,A maximal list of indices i such that the set of all sectionlist[i] is independent,0,2,0,0,0,0,0,0,0,AdelMat,,0,0,82,,82,-38,-38,-38,-38,-38
S,ExtractCoeffs,"Given a matrix mat defined over a Laurent series ring with uniformizer pi, return the matrix of the coefficients of pi^i in each component",0,2,0,0,0,0,0,0,0,148,,0,0,-1,,177,-38,-38,-38,-38,-38
S,TwistAndCompare,Operates on admat1 and admat2 in place,0,3,0,0,0,0,0,0,0,148,,0,0,AdelMat,,0,0,AdelMat,,36,AdelMat,AdelMat,-38,-38,-38
S,TwistForEnoughSections,Test,0,2,0,0,0,0,0,0,0,AdelMat,,0,0,AdelMat,,36,AdelMat,AdelMat,82,82,82
S,IsIsomorphic,"return true if admat1 and admat2 define isomorphic vector bundles, false if not. If true, also return the divisor D and matrix m in GL2(F) such that admat1 and m * admat2 * D are equivalent mod GL2(O) (ie up to ""Reduce"")",0,2,0,0,0,0,0,0,0,AdelMat,,0,0,AdelMat,,36,82,-38,-38,-38,-38
