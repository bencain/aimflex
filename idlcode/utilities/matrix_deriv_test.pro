pro matrix_deriv_test

n=2;4

xs=[0d,1d];,1d,0d]
ys=[0d,1d];,1d,1d]

pos=[xs,ys]

shear=dblarr(2*n)

shear[0:n-1]=3d*pos[0:n-1]+4d*pos[n:2*n-1]
shear[n:2*n-1]=-pos[0:n-1]-2d*pos[n:2*n-1]

Dx=dblarr(2*n,2*n)
Dy=dblarr(2*n,2*n)

for i=0,2*n-1 do begin
   for j=0,2*n-1 do begin
      
      Dx[i mod n,j mod n]=1d/(pos[i mod n]-pos[j mod n])
      Dy[i mod n,j mod n]=1d/(pos[n+(i mod n)]-pos[n+(j mod n)])
      Dx[n+(i mod n),n+(j mod n)]=1d/(pos[i mod n]-pos[j mod n])
      Dy[n+(i mod n),n+(j mod n)]=1d/(pos[n+(i mod n)]-pos[n+(j mod n)])

   endfor
endfor

print,Dx
print,''
print,Dy

end
