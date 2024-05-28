function F_bil = get_Fbil(x,y,c,x0,y0)

if c > 0
    t0 = y0/x0;
    F_bil = 0.5*(t0*x^2 + 1/t0*y^2);
else
    F_bil = 0.25*(x-y)^2 + 0.25*(x0+y0)^2 - 0.5*(x0+y0)*(x+y);
end

end % EOF