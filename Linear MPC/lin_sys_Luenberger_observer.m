function  [y1,x1_updated] = lin_sys_Luenberger_observer(x0,sys,u0,y_plant_0,K1)
x1_predicted    = sys.a*x0 + sys.b*u0;
y_0             = sys.c*x0;
e_0             = y_plant_0 - y_0;
x1_updated      = x1_predicted + K1*e_0;
y1              = sys.c*x1_updated;
end