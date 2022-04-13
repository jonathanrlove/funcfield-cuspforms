174,0
S,SmallerSlant,Return an equivalent adelic matrix with slant degree at most 2g,0,1,0,0,0,0,0,0,0,AdelMat,,AdelMat,-38,-38,-38,-38,-38
S,HeckeNeighbors,return the neighbors of adelmat according to the Hecke operator at v,0,2,0,0,0,0,0,0,0,232,,0,0,AdelMat,,82,-38,-38,-38,-38,-38
S,RatPlaces,return a list of rational places of C,0,1,0,0,0,0,0,0,0,371,,82,-38,-38,-38,-38,-38
S,HigherDegreePlace,Returns a place of minimal degree greater than 1,0,1,0,0,0,0,0,0,0,371,,232,-38,-38,-38,-38,-38
T,HeckeInfo,-,0
A,HeckeInfo,10,curve,alldoubcos,allseenreps,slants,terminals,heckes,heckeplaces,spaces,eigenforms,eisensteins
S,NewHeckeInfo,Creates a new HeckeInfo,0,1,0,0,0,0,0,0,0,371,,HeckeInfo,-38,-38,-38,-38,-38
S,HeckeSearch,"Hecke search. To do: predict number of rows/columns needed rather than just using an absurd number. double cosets with slant less than 2-2g have the property that any cuspidal function vanishes on them; to compute cuspidal functions, we do not need to record the neighbors of such double cosets that are themselves of smaller slant. Optional parameter terminalextension computes all neighbors for double cosets for slant at least 2-2g-terminalextension. status can be 0, 1, or 2, depending on how many intermediate results to be printed as the computation progresses",0,2,0,0,1,0,0,0,0,82,,0,0,HeckeInfo,,-38,-38,-38,-38,-38,-38
S,HeckeMatrix,Print the Hecke matrix at a given place,0,2,0,0,0,0,0,0,0,232,,0,0,HeckeInfo,,177,-38,-38,-38,-38,-38
S,SimulEigens,"Computes simultaneous eigenvectors for a sequence of matrices acting on a space. If space can be decomposed into simultaneous eigenspaces, return true and a sequence consisting of vectors and eigenvalues. Otherwise return false",0,2,0,0,0,0,0,0,0,82,,0,0,191,,168,-38,-38,-38,-38,-38
S,ComputeEigenforms,"Compute eigenforms of heckeinfo, and assign them to heckeinfo`eigenforms",0,1,0,0,1,0,0,0,0,HeckeInfo,,-38,-38,-38,-38,-38,-38
S,IsSimultaneousEigenform,Returns whether a form is an eigenform for Hecke operators at all places where heckeinfo has been computed,0,2,0,0,0,0,0,0,0,HeckeInfo,,0,0,160,,36,82,82,-38,-38,-38
S,EigenformData,Returns whether a form is an eigenform for Hecke operators at all places where heckeinfo has been computed,0,1,0,0,0,0,0,0,0,HeckeInfo,,82,-38,-38,-38,-38,-38
