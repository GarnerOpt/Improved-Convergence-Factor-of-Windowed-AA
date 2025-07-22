%==========================================================================
% This script compares Anderson Acceleration with full-memory, with
% restarting, and with a sliding frame. 
%==========================================================================
clc;
clear all;
close all;

% Dimension
n = 500;
maxiter = 1000; 
tol = 1e-12; 
beta = 1;

% Data generation and setting random number generator
rng(13)
A = randn(n,n);
A = A + A';
A = A/(1.01*norm(A,2));
%[A,~] = linear_operator_data(n,0,1,-1);

xstar = rand(n,1);
x0 = rand(n,1);
b = A*xstar;
q_data = {A,b};

%--------------------------------------------------------------------------
% Full-Memory
%--------------------------------------------------------------------------
m = maxiter;
[xfinal_full, x_iter_full, err_iter_full, runtime_full] = AA_Rn(@q, q_data, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_full)
    full_err(i) = norm(x_iter_full{i}-(A*x_iter_full{i}+b),2);
end

%--------------------------------------------------------------------------
% Sliding Window of length m
%--------------------------------------------------------------------------
m = 2;  
[xfinal_slide, x_iter_slide, err_iter_slide, runtime_slide] = AA_Rn(@q, q_data, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_slide)
    slide_err(i) = norm(x_iter_slide{i}-(A*x_iter_slide{i}+b),2);
end

%--------------------------------------------------------------------------
% With Restarting 
%--------------------------------------------------------------------------
m = 2;
[xfinal_restart, x_iter_restart, err_iter_restart, runtime_restart] = AArestart_Rn(@q, q_data, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_restart)
    restart_err(i) = norm(x_iter_restart{i}-(A*x_iter_restart{i}+b),2);
end

%--------------------------------------------------------------------------
% Plot Error of Each Method
%--------------------------------------------------------------------------
figure(1);
semilogy(full_err,'r-','LineWidth',1.5);
hold on;
semilogy(slide_err,'b-','LineWidth',1.5);
hold on;
semilogy(restart_err,'k-','LineWidth',1.5);
title('Comparing AA Variants (Linear)','FontSize',14,'FontWeight','bold');
ylabel('Fixed-Point Error','FontSize',14,'FontWeight','bold');
xlabel('Iterations','FontSize',14,'FontWeight','bold');
legend({'Full-Memory','AA(m)-Sliding','AA(m)-Restart'},'FontSize',14,'FontWeight','bold');


%%%
% Helper Functions
%%%

function out = q(x,q_data)
    out = q_data{1}*x + q_data{2};
end
