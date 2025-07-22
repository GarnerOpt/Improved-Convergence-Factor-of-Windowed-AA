%==========================================================================
% AISTATS TME Experiment 4
%==========================================================================

%--------------------------------------------------------------------------
% Parameters for Algorithms 
%--------------------------------------------------------------------------
max_iter = 5000;
tol = 1e-12;
beta = 1;

%--------------------------------------------------------------------------
% Data Generation
%--------------------------------------------------------------------------
total_tests = 10;
total_models = 10;
p = 200;
n = 210;
type = 1;

% ---------------------------------------------------------------------
figure;
maxAA_iters=0;
Iterations = [];
Methods = [];
Times = [];

tests = 0;
for data_ins = 1:total_models
    rng(data_ins);
    [X,p,n] = TME_data(p,n,type);

    for test_no = 1:total_tests
        tests = tests + 1;
        tests
        %======================================
        % Fixed-Point Approach
        %======================================
        temp = randn(p,n);
        cov_0=(p/n)*temp*temp';
        [cov_sol,cov_iter,temp,err_iter_FP,iter,time_FP] = m_estimator(X',cov_0,max_iter,tol);

        r_iter_fp = r_linear_err(cov_sol,cov_iter);
        FP_iters(tests) = length(err_iter_FP);
        R_rate_FP{tests} = r_iter_fp;
        Err_FP{tests} = err_iter_FP;
    
        Methods(end+1)=1;
        Iterations(end+1)=FP_iters(end);
        Times(end+1)=time_FP;
        
        %======================================
        % AA(m) on Symmetric TME Formulation
        %====================================== 
        w_0 = randn(n,1);
        m_values = [1,2,3];
        for l=1:length(m_values)
            m = m_values(l);
            [w_final, w_iter, AA_err_iter, runtime_AA] = AA_Rn(@q, X, w_0, m, beta, max_iter, tol);

            Methods(end+1)= m+1;
            Iterations(end+1)=length(AA_err_iter);
            Times(end+1)=runtime_AA;
    
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
                plot([1,length(err_iter_FP)],[RHS,RHS],'r-','LineWidth',1.5);
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
end

%==========================================================================
% Make Box-Plot of Numerial Experiments
%==========================================================================
indx=find(Methods==1);
Times_AA = Times;
Times_AA(indx)=[];

Iterations_AA = Iterations;
Iterations_AA(indx)=[];

Methods_AA = Methods;
Methods_AA(indx)=[];

Alg_AA = categorical(Methods_AA,[2,3,4],{'AA(1)','AA(2)','AA(3)'});
Alg = categorical(Methods,[1,2,3,4],{'FP','AA(1)','AA(2)','AA(3)'});

%--------------------------------------------------------------------------
% Boxplot for Compute Time
%--------------------------------------------------------------------------
figure;
b1 = axes();
boxplot(b1,Times,Alg)
title('Compute Time Comparison','FontSize',16,'FontWeight','bold')
xlabel('Method','FontSize',16,'FontWeight','bold')
ylabel('Time (s)','FontSize',16,'FontWeight','bold')
ax3 = gca;
ax3.FontWeight = 'bold';

%--------------------------------------------------------------------------
% Boxplot for Total Iterations
%--------------------------------------------------------------------------
figure;
b1 = axes();
boxplot(b1,Iterations,Alg)
title('Total Iteration Comparison','FontSize',16,'FontWeight','bold')
xlabel('Method','FontSize',16,'FontWeight','bold')
ylabel('Iterations','FontSize',16,'FontWeight','bold')
ax1 = gca;
ax1.FontWeight = 'bold';

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
    fprintf('\n_________________________________________\n')
    fprintf('_________________________________________\n')
    fprintf('Standard TME Fixed-Point Iteration: \n')
    fprintf('\t Final Error Value: %7.4e\n',err_iter(end))
    fprintf('\t Total Iterations: %5.0f\n',iter-1)
    fprintf('\t Total time: %5.2f seconds \n',time);
end

function wnew = q(w,Xdata)
    [p,n]=size(Xdata);
    wnew = zeros(n,1);
    delta = eps;
    M    = inv(Xdata*diag(exp(w))*Xdata'+ delta*eye(p));
    wnew = -log((n/p)*diag(Xdata'*M*Xdata));
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
