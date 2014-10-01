function rotate_image, image, angle

dims=size(image,/dimensions)

ctr=double((dims-1)/2)

; We want (x' - x_c) = R(phi)*(x - x_c)

p=dblarr(2,2)
p[0,0]=(1d - cos(angle))*ctr[0] - sin(angle)*ctr[1]
p[1,1]=0d
p[1,0]=sin(angle)
p[0,1]=cos(angle)

q=dblarr(2,2)
q[0,0]=sin(angle)*ctr[0] + (1d - cos(angle))*ctr[1]
q[1,1]=0d
q[1,0]=cos(angle)
q[0,1]=-sin(angle)

rotated_image=poly_2d(image,p,q,missing=0d)

return,rotated_image

end
