function [A_eq,Q_X,f] = constant_term(sys)
A_eq = zeros((sys.n_x+sys.n_u+sys.n_y)*sys.h,(sys.n_x+2*sys.n_u+sys.n_y)*sys.h);
A_eq(1:sys.n_x+sys.n_u,1:sys.n_x+2*sys.n_u) = [-sys.bp,eye(sys.n_x+sys.n_u)];
Q_X_1 = blkdiag(sys.R_du,sys.Q_xp);
Q_X_2 = sys.Q_y;
A_eq((sys.n_x+sys.n_u)*sys.h+1:(sys.n_x+sys.n_u)*sys.h+sys.n_y,sys.n_u+1:sys.n_u+(sys.n_x+sys.n_u)) = -sys.cp;
for i=2:sys.h
    A_eq((i-1)*(sys.n_x+sys.n_u)+1:(i)*(sys.n_x+sys.n_u),sys.n_u+(i-2)*(sys.n_x+2*sys.n_u)+1:sys.n_u+(i-1)*(sys.n_x+2*sys.n_u)+sys.n_x+sys.n_u) = [-sys.ap,-sys.bp,eye(sys.n_x+sys.n_u)];    
    Q_X_1 = blkdiag(Q_X_1,sys.R_du,sys.Q_xp);
    Q_X_2 = blkdiag(Q_X_2,sys.Q_y);
    A_eq((sys.n_x+sys.n_u)*sys.h+(i-1)*sys.n_y+1:(sys.n_x+sys.n_u)*sys.h+i*sys.n_y,sys.n_u+(i-1)*(sys.n_x+2*sys.n_u)+1:sys.n_u+(i-1)*(sys.n_x+2*sys.n_u)+sys.n_x+sys.n_u) = -sys.cp;
end
Q_X = blkdiag(Q_X_1,Q_X_2);
A_eq((sys.n_x+sys.n_u)*sys.h+1:end,(sys.n_x+2*sys.n_u)*sys.h+1:end) = eye(sys.h*sys.n_y);
f = zeros(size(Q_X,1),1);
end