%==========================================================================
% AISTATS TME Experiment 1
%==========================================================================
% GOAL: Demonstrate the accuracy of the theory for AA(m) applied in the
% nonlinear operator setting with the TME. 
%==========================================================================

fprintf('\n=============================================================\n')
fprintf('=============================================================\n')
fprintf('               NONLINEAR OPERATOR EXPERIMENT 1                 \n')
fprintf('=============================================================\n')
fprintf('=============================================================\n\n')

%--------------------------------------------------------------------------
% Parameters for Algorithms 
%--------------------------------------------------------------------------
max_iter = 4000;
tol = 1e-12;
beta = 1;

%--------------------------------------------------------------------------
% Data Generation and fixed randome number generator
%--------------------------------------------------------------------------
p = 100;
n = 110;
type = 1;
rng(3)
[X,p,n] = TME_data(p,n,type);

% Exact solution
[cov_sol,cov_iter,temp,err_iter,iter,time_FP] = m_estimator(X',(p/n)*(X*X'),max_iter,tol);

total_tests = 100;

maxAA_iters=0;
for tests = 1:total_tests

    w0 = randn(n,1);
    %======================================
    % TME Alt. Fixed-Pt 
    %======================================

    [w_sol_FP, w_iter_FP, err_iter_FP, time_FP] = TME_FP_Sym(X,w0,max_iter,tol);

    % r-linear convergence estimation
    w_new_FP = {};

    % transform each vector w into scaled matrix
    for i=1:length(w_iter_FP)
        w = w_iter_FP{i};
        Qw = (p/n)*X*diag(exp(w))*X';
        Qw = p/trace(Qw)*Qw; 
        w_new_FP{i} = Qw;
    end

    % compute estimate
    r_iter_FP = r_linear_err(cov_sol,w_new_FP);

    % store FP information
    R_rate_FP{tests} = r_iter_FP;
    Err_FP{tests} = err_iter_FP;
    FP_iters(tests) = length(err_iter_FP);
    FP_sol_diff(tests) = norm(cov_sol-w_new_FP{end});

    %======================================
    % AA(m) on Symmetric TME Formulation
    %====================================== 

    m_values = [1,2,3];
    for l=1:length(m_values)
        m = m_values(l);
        [w_final, w_iter, AA_err_iter, runtime_AA] = AA_Rn(@q, X, w0, m, beta, max_iter, tol);

        % Compute the information for the est. of r-linear converg. rate
        if (tests == 1) && (m==1)
            W = q_grad(w_final,X);
            eigW = sort(eig(W)); 
            maxW = eigW(end-1);
            minW = eigW(1);
            normW = maxW;
            aval = linspace(0,1e3,10e7);
            for i=1:length(aval)
                a = aval(i);
                y(i) = r_bound(a,maxW,minW);
            end
            theta = max(y);
            r_est = sqrt(normW*(sin(theta)));
            RHS = sin(theta)*normW;
            %plot([1,length(err_iter)],[RHS,RHS],'r-','LineWidth',1.5);
        end
    
        w_new = {};
        for i=1:length(w_iter)
            w = w_iter{i};
            Qw = (p/n)*X*diag(exp(w))*X';
            Qw = p/trace(Qw)*Qw; 
            w_new{i} = Qw;
        end
    
        
        for i=2:length(w_iter)-1
            temp_n = norm(w_iter{i+1} - q(w_iter{i+1},X));
            temp_d = norm(w_iter{i-1} - q(w_iter{i-1},X));
            lhs(i-1) = temp_n/temp_d;
        end
        LHS{l,tests} = lhs;
        lhs=[];
    
        r_iter_AA = r_linear_err(cov_sol,w_new);
        AA_sol_diff{l,tests}=norm(cov_sol-w_new{end});
        R_rate_AA{l,tests} = r_iter_AA;
        Err_AA{l,tests} = AA_err_iter;
        if length(w_iter)>maxAA_iters
            maxAA_iters = length(w_iter);
        end
    end
end

%======================================
% R-Linear Rate Estimate Plot
%======================================
% Plot Results: R-Linear Convergence Rate
figure;
for i=1:total_tests
    plot(R_rate_FP{i},'b-','LineWidth',1.5);
    hold on;
    plot(R_rate_AA{1,i},'r-','LineWidth',1.5);
    hold on;
    plot(R_rate_AA{2,i},'k-','LineWidth',1.5);
    hold on;
    plot(R_rate_AA{3,i},'g-','LineWidth',1.5);
    hold on;
    plot([1,length(Err_FP{1})],[r_est,r_est],'k--','LineWidth',1.5);
    hold on;
end
legend({'FP','AA(1)','AA(2)','AA(3)','Thm 2 Est.'},'FontSize',14,'FontWeight','bold','Location','best');
set(gca, 'YDir', 'reverse')
xlabel('Iterations (k)','FontSize',16,'FontWeight','bold');
ylabel('r_{est}(k)','FontSize',16,'FontWeight','bold')
axis([0,maxAA_iters+10,0.55, 1.25])
title('Comparing R-Linear Convergence Rates','FontSize',16,'FontWeight','bold');


%======================================
% R-Linear Rate Ratio Estimate
%======================================
figure;
for i=1:total_tests
    plot([1,length(Err_FP{1})],[RHS,RHS],'k--','LineWidth',2);
    hold on;
    plot(LHS{1,i},'r-','LineWidth',1.5);
    hold on;
    plot(LHS{2,i},'b-','LineWidth',1.5);
    hold on;
    plot(LHS{3,i},'g-','LineWidth',1.5);
    hold on;
end
legend({'RHS','LHS-AA(1)','LHS-AA(2)','LHS-AA(3)'},'FontSize',14,'FontWeight','bold','Location','best');
xlabel('Iterations (k)','FontSize',16,'FontWeight','bold');
ylabel('Values in Eq. (17)','FontSize',16,'FontWeight','bold')
axis([0,maxAA_iters+10,0, 1.05*RHS])
title('Tightness of Ratio','FontSize',16,'FontWeight','bold');


%--------------------------------------------------------------------------
% Helper Functions
%--------------------------------------------------------------------------

function r_iter = r_linear_err(Xsol,Xiter)
    r_iter = zeros(length(Xiter),1);
    for i=1:length(r_iter)
        r_iter(i) = norm(Xsol-Xiter{i},2)^(1/i);
    end
    r_iter = max_n_to_end(r_iter);
end

function r_err = r_linear_err_vec(x_sol, x_iter)
    for i=1:length(x_iter)
        r_err(i) = norm(x_sol - x_iter{i})^(1/i);
    end
end

function new_vec = max_n_to_end(vec)
    new_vec = zeros(length(vec),1);
    for i=1:length(vec)
        new_vec(i) = max(vec(i:end));
    end
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
%     fprintf('\t Final Error Value: %7.4e\n',err_iter(end))
%     fprintf('\t Total Iterations: %5.0f\n',iter-1)
%     fprintf('\t Total time: %5.2f seconds \n',time);
end

function wnew = q(w,Xdata)
    [p,n]=size(Xdata);
    wnew = zeros(n,1);
    delta = eps;
    M    = inv(Xdata*diag(exp(w))*Xdata'+ delta*eye(p));
    wnew = -log((n/p)*diag(Xdata'*M*Xdata));
end

function [w_sol,w_iter,err_iter,time]=TME_FP_Sym(X,w0,max_iter,tol)
    k=1;
    tic;
    wnew = q(w0,X);
    err_iter(1)=norm(w0-wnew,2);
    w_iter{1}= wnew;
    w0=wnew;
    while (err_iter(end)>tol)&&(k<=max_iter)
        wnew = q(w0,X);
        w_iter{k+1}=wnew;
        err_iter(end+1) = norm(w0-wnew,2);
        w0 = wnew;
        k = k+1;
    end
    w_sol = wnew;
    time = toc;
    fprintf('\n_________________________________________\n')
    fprintf('_________________________________________\n')
    fprintf('Alternative TME Fixed-Point Iteration: \n')
    fprintf('\t Final Error Value: %7.4e\n',err_iter(end))
    fprintf('\t Total Iterations: %5.0f\n',k-1)
    fprintf('\t Total time: %5.2f seconds \n',time);
end

function W = q_grad(w,Xdata)
    [p,n]=size(Xdata);
    W = zeros(n,n);
    delta = eps;
    Sigma = inv((p/n)*(Xdata*diag(exp(w))*Xdata'+ delta*eye(p)));
    for j=1:n
        xj = Xdata(:,j);
        for i=j:n
            xi = Xdata(:,i);
            W(j,i) = exp(w(j))*xj'*Sigma*((p/n)*exp(w(i))*(xi*xi'))*Sigma*xj;
            W(i,j) = W(j,i);
        end
    end
end

function out = r_bound(a,maxW,minW)
    num = abs(maxW-minW)*a;
    den = sqrt(4+((2-maxW-minW)^2)*(a^2));
    term1 = asin(num/den);
    term2 = abs( atan(a) - atan((1-0.5*(maxW + minW))*a)  );
    out = term1 + term2;  
end

