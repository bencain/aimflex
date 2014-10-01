function convert_sym_matrix, matrix

dims=size(matrix,/dimensions)

; Check if we have a linear input - then output a 2D version
if (n_elements(dims) eq 1) or (dims[0] eq 1) then begin
   n1D=n_elements(matrix)
   n2D=long((1d +sqrt(1d + 8d*n1D))/2d)
   out_matrix=dblarr(n2D,n2D)

   k=0
   for j=0,n2D-2 do for i=j+1,n2D-1 do begin
;      k=i-1+(n2D-1)*j
      out_matrix[i,j]=matrix[k]
      out_matrix[j,i]=matrix[k]
      k++
   endfor
   for i=0,n2D-1 do out_matrix[i,i]=1d

; Otherwise, shift a 2D into a 1D output.
endif else begin
   n2D=min(dims[0:1])
   n1D=(n2D^2-n2D)/2
   out_matrix=dblarr(n1D)

   k=0
   for j=0,n2D-2 do for i=j+1,n2D-1 do begin
;      k=i-1+(n2D-1)*j
      out_matrix[k]=matrix[i,j]
      k++
   endfor
endelse

return,out_matrix
end
