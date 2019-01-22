function b_eq = b_eq_calc(sys,x_spp,x_0,u__1,y_sp)
b_eq = zeros(sys.h*(sys.n_x+sys.n_u+sys.n_y),1);
b_eq(1:sys.n_x+sys.n_u) = sys.ap*[x_0;u__1]-x_spp(:,1);
i = 1;
b_eq(sys.h*(sys.n_x+sys.n_u)+(i-1)*sys.n_y+1:sys.h*(sys.n_x+sys.n_u)+i*sys.n_y) = +sys.cp*x_spp(:,i)-y_sp(:,i);
for i=2:sys.h
    b_eq((i-1)*(sys.n_x+sys.n_u)+1:i*(sys.n_x+sys.n_u)) = sys.ap*x_spp(:,i-1)-x_spp(:,i);
    b_eq(sys.h*(sys.n_x+sys.n_u)+(i-1)*sys.n_y+1:sys.h*(sys.n_x+sys.n_u)+i*sys.n_y) = +sys.cp*x_spp(:,i)-y_sp(:,i);
end
end