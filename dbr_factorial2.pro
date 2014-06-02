;dbr_factorial2.pro

function dbr_factorial2,n

start = n mod 2
if start eq 0 then start=2

num = start
factorial2 = start
while num le n-2 do begin
   num = num+2
   factorial2 *= num
endwhile 
return,factorial2
end
