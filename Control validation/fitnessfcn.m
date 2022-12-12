function J = fitnessfcn(X, n,A,B,C,D,settings,r0,v0,A_BN0,w0)


Kp = X(1:3);
Kd = X(4:6);
test = sim("PE_test.slx",'SrcWorkspace','Current');
w = reshape(test.w', [1, size(test.w, 1)*size(test.w, 2)]);
w_ref = zeros(1, length(w));
w_ref(2:3:end) = n;
e = reshape(test.e', [1, size(test.e, 1)*size(test.e, 2)]);
u = reshape(test.u', [1, size(test.u, 1)*size(test.u, 2)]);

J = 1/2  * (w-w_ref)*(w-w_ref)' + 1/2 * e*e' * 1/2 * u*u' ;