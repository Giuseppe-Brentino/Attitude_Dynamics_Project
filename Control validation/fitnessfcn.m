function J = fitnessfcn(X,settings,sat,environment,n,r0,v0,A0,A_BN0)


K_cx = X(1:3);
K_cz = X(4:6);
K_w = X(7);
test = sim("POINTING_2.slx",'SrcWorkspace','Current');

J = trapz(test.ew + test.theta_x + test.theta_z);