function J = fitnessfcn(X, n)



test = sim("PE_test.slx");
w = reshape(test.w', [1, size(test.w, 1)*size(test.w, 2)]);
w_ref = zeros(1, length(w));
w_ref(2:3:end) = n;
q = reshape(test.q_e', [1, size(test.q_e, 1)*size(test.q_e, 2)]);
Mc = reshape(test.Mc', [1, size(test.Mc, 1)*size(test.Mc, 2)]);

J = 1/2  * ( )^2;