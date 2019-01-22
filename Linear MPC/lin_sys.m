function  [y,x] = lin_sys(x0,sys,u)
x = sys.a*x0 + sys.b*u;
y = sys.c*x + sys.d*u;
end