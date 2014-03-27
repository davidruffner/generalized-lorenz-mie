;script for making a cross section

function fieldcrossectionpoints, r0,dv0,dv1,n0,n1

pts = fltarr(3,n0,n1)

for i=0,n0-1 do begin
   for j = 0,n1-1 do begin
      point = r0+dv0*double(i)+dv1*double(j)
      pts[*,i,j] = point
   endfor
endfor

return, pts
end

