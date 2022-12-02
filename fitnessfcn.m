function J = fitnessfcn(X,A,B,final_w,high_torque,Ts,sat,settings)

A_tot = A + X(1)*eye(3);
Q = diag(X(2)*ones(3,1));
R = diag(X(3)*ones(3,1));
k = lqr(A_tot,B,Q,R);
test = sim("controlTuning.slx",'SrcWorkspace','Current');


max_torque = max(max(abs(test.torque)));

J = 1/2  * 5 * ( (test.w(end,:)-final_w')*(test.w(end,:)'-final_w) ) + 1/2 * 1 * (max_torque-high_torque)^2;