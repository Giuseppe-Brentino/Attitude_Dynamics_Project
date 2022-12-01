function J = fitnessfcn(X,I,w,sat,settings)

A =[0 (I(2)-I(3))/I(1)*w (I(2)-I(3))/I(1)*w;...
    (I(3)-I(1))/I(2)*w 0 (I(3)-I(1))/I(2)*w;...
    (I(1)-I(2))/I(3)*w (I(1)-I(2))/I(3)*w 0] + X(1)*eye(3);
B = inv(diag(I));
Q = diag(X(2)*ones(3,1));
R = diag(X(3)*ones(3,1));
k = lqr(A,B,Q,R);
test = sim("controlTuning.slx",'SrcWorkspace','Current');

for i = 1:size(test.torque,1)
    torque(i) = norm(test.torque(i,:));
end
max_torque = max(torque);
J = test.w(end,:)*test.w(end,:)' + max_torque;