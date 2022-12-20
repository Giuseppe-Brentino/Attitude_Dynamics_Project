function J = fitnessfcn(X,settings,sat,environment,thruster,n,r0,v0,A0,A_BN0)


thrusters.K_m = X(1);
thrusters.K_m = X(2);

test = sim("PE_test.slx",'SrcWorkspace','Current');

if isfield('test', 'ae')
    J = trapz(test.ae(:,1) + test.ae(:,2) + test.ae(:,3));
else
    J=10^10;
end