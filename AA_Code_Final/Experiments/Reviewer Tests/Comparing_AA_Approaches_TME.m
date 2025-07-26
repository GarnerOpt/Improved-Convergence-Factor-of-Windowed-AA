%==========================================================================
% This script compares Anderson Acceleration with full-memory, with
% restarting, and with a sliding frame on the two TME Data Models. 
%==========================================================================
clc;
clear all;
close all;

%==========================================================================
% TME Data Model 1
%==========================================================================
%--------------------------------------------------------------------------
% Algorithm Settings and Data Itialization
%--------------------------------------------------------------------------
 rng(5);
 type = 1;
 p = 100;
 n = 110;
 [X,p,n] = TME_data(p,n,type);
 x0 = randn(n,1);
 maxiter = 5000;
 tol = 10^(-12);
 beta = 1;

%--------------------------------------------------------------------------
% Baseline Solution w/ STD Fixed-Point Approach
%--------------------------------------------------------------------------
temp = randn(p,n);
cov_0=(p/n)*temp*temp';
[cov_sol,cov_iter,temp,err_iter_FP,iter,time_FP] = m_estimator(X',cov_0,maxiter,tol);

%--------------------------------------------------------------------------
% Full-Memory
%--------------------------------------------------------------------------
m = maxiter;
[xfinal_full, x_iter_full, err_iter_full, runtime_full] = AA_Rn(@q, X, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_full)
    w_full = Form_TME_Sol(x_iter_full{i},X,p,n);
    full_matrix_err(i) = norm(w_full-cov_sol);
end
W_full = Form_TME_Sol(xfinal_full,X,p,n);

%--------------------------------------------------------------------------
% Sliding Window of length m
%--------------------------------------------------------------------------
m = 2;
[xfinal_slide, x_iter_slide, err_iter_slide, runtime_slide] = AA_Rn(@q, X, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_slide)
    w_slide = Form_TME_Sol(x_iter_slide{i},X,p,n);
    slide_matrix_err(i) = norm(w_slide-cov_sol);
end

%--------------------------------------------------------------------------
% With Restarting 
%--------------------------------------------------------------------------
m = 2; 
[xfinal_restart, x_iter_restart, err_iter_restart, runtime_restart] = AArestart_Rn(@q, X, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_restart)
    w_restart = Form_TME_Sol(x_iter_restart{i},X,p,n);
    restart_matrix_err(i) = norm(w_restart-cov_sol);
end

%--------------------------------------------------------------------------
% Plot Error of Each Method
%--------------------------------------------------------------------------
figure;
semilogy(err_iter_full,'r-','LineWidth',1.75);
hold on;
semilogy(err_iter_slide,'b-','LineWidth',1.75);
hold on;
semilogy(err_iter_restart,'k-','LineWidth',1.75);
title('TME Data Model 1','FontSize',14,'FontWeight','bold');
ylabel('Fixed-Point Error','FontSize',14,'FontWeight','bold');
xlabel('Iterations','FontSize',14,'FontWeight','bold');
legend({'AA(\infty)','AA(2)','restart-AA(2)'},'FontSize',14,'FontWeight','bold');

%==========================================================================
%==========================================================================
% TME Data Model 2 
%==========================================================================

%--------------------------------------------------------------------------
% Algorithm Settings and Data Itialization
%--------------------------------------------------------------------------
 rng(5);
 [X,p,n] = TME_data(0,0,2);
 x0 = randn(n,1);
 maxiter = 5000;
 tol = 10^(-12);
 beta = 1;

%--------------------------------------------------------------------------
% Baseline Solution w/ STD Fixed-Point Approach
%--------------------------------------------------------------------------
temp = randn(p,n);
cov_0=(p/n)*temp*temp';
[cov_sol,cov_iter,temp,err_iter_FP,iter,time_FP] = m_estimator(X',cov_0,maxiter,tol);

%--------------------------------------------------------------------------
% Full-Memory
%--------------------------------------------------------------------------
m = maxiter;
[xfinal_full, x_iter_full, err_iter_full, runtime_full] = AA_Rn(@q, X, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_full)
    w_full = Form_TME_Sol(x_iter_full{i},X,p,n);
    full_matrix_err(i) = norm(w_full-cov_sol);
end
W_full = Form_TME_Sol(xfinal_full,X,p,n);

%--------------------------------------------------------------------------
% Sliding Window of length m
%--------------------------------------------------------------------------
m = 2;
[xfinal_slide, x_iter_slide, err_iter_slide, runtime_slide] = AA_Rn(@q, X, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_slide)
    w_slide = Form_TME_Sol(x_iter_slide{i},X,p,n);
    slide_matrix_err(i) = norm(w_slide-cov_sol);
end

%--------------------------------------------------------------------------
% With Restarting 
%--------------------------------------------------------------------------
m = 2; 
[xfinal_restart, x_iter_restart, err_iter_restart, runtime_restart] = AArestart_Rn(@q, X, x0, m, beta, maxiter, tol);
for i=1:length(x_iter_restart)
    w_restart = Form_TME_Sol(x_iter_restart{i},X,p,n);
    restart_matrix_err(i) = norm(w_restart-cov_sol);
end

%--------------------------------------------------------------------------
% Plot Error of Each Method
%--------------------------------------------------------------------------
figure;
semilogy(err_iter_full,'r-','LineWidth',1.75);
hold on;
semilogy(err_iter_slide,'b-','LineWidth',1.75);
hold on;
semilogy(err_iter_restart,'k-','LineWidth',1.75);
title('TME Data Model 2','FontSize',14,'FontWeight','bold');
ylabel('Fixed-Point Error','FontSize',14,'FontWeight','bold');
xlabel('Iterations','FontSize',14,'FontWeight','bold');
legend({'AA(\infty)','AA(2)','restart-AA(2)'},'FontSize',14,'FontWeight','bold');

%==========================================================================
%%%
% Helper Functions
%%%
%==========================================================================

function wnew = q(w,Xdata)
    [p,n]=size(Xdata);
    wnew = zeros(n,1);
    delta = eps;
    M    = inv(Xdata*diag(exp(w))*Xdata'+ delta*eye(p));
    wnew = -log((n/p)*diag(Xdata'*M*Xdata));
end

function w_sol = Form_TME_Sol(w,X,p,n)
        Qw = (p/n)*X*diag(exp(w))*X';
        Qw = p/trace(Qw)*Qw; 
        w_sol = Qw;
end
  
function [cov_sol,cov_iter,temp,err_iter,iter,time]=m_estimator(X,cov_0,max_iter,tol)
    % X: data matrix with each row representing a point
    % cov: the estimated covariance
    tic;
    [N,D]=size(X);
    initcov=eye(D);
    oldcov=initcov-1;
    cov=cov_0;
    cov_iter{1}=cov;
    err=1;
    iter=1;
    delta=eps;%regulirization parameter
    while (err >tol) && (iter< max_iter)
        temp=X*(cov+delta*eye(D))^-1;   
        d=sum(temp.*conj(X),2);
        oldcov=cov;
        temp=((real(d)+delta*ones(N,1)).^-1);
        cov=X'.*repmat(((temp)),1,D)'*X/N;
        cov=D*(cov/trace(cov));
        err = norm(oldcov-cov,'fro');
        err_iter(iter)=err;
        iter=iter+1;
        cov_iter{iter}=cov;
    end
    cov_sol = cov;
    time=toc;
    fprintf('\n_________________________________________\n')
    fprintf('_________________________________________\n')
    fprintf('Standard TME Fixed-Point Iteration: \n')
    fprintf('\t Final Error Value: %7.4e\n',err_iter(end))
    fprintf('\t Total Iterations: %5.0f\n',iter-1)
    fprintf('\t Total time: %5.2f seconds \n',time);
end