function F_qua = get_Fqua(x,c,x0)

F_qua = 2*(c - x0)'*(x - x0) - norm(x0 - c)^2;

end % EOF