function eigenvals, M, evecs=evecs
  
; This function returns eigenvectors and eigenvalues from the singular
; value decomposition of the matrix.  The eigenvalues are returned in
; an N element array EVALS where M is an NxN matrix. The eigenvectors
; are returned in an NxN array EVECS where the vector EVECS[i,*]
; corresponds to the eigenvalue EVALS[i].

svdc,M,evals,U,evecs,/double

return,evals

end
