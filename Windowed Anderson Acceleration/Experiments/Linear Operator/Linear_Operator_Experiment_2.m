%==========================================================================
% Linear Operator - Verying Eq. (4) Tightness Test 2 
%==========================================================================
% GOAL: 
% Empirically verying the tightness of Eq.(4) for the case M is 5x5
% non-diagonal, symmetric. 
%--------------------------------------------------------------------------

fprintf('\n=============================================================\n')
fprintf('=============================================================\n')
fprintf('               LINEAR OPERATOR EXPERIMENT 2                   \n')
fprintf('=============================================================\n')
fprintf('=============================================================\n\n')

% Dimension
n = 500;
maxiter = 1000; 
tol = 1e-12; 
m = 2;  
beta = 1;

% Data generation and setting random number generator
rng(13)
[M,~] = linear_operator_data(n,0,1,-1);

% Computing Constants for the estimate from Theorem 2 
W = beta*M + (1-beta)*eye(n,n);
eigW = eig(W); 
maxW = max(eigW);
minW = min(eigW);
aval = linspace(0,1e3,10e7);
for i=1:length(aval)
    a = aval(i);
    y(i) = r_bound(a,maxW,minW);
end
theta = max(y);
r_est = sqrt(norm(W,2)*abs(sin(theta)));


% Plot Figure to Show Estimate from Theorem 2 is tight
figure;
rhs = norm(W,2)*sin(theta);
min_diff = inf;
plot([1,100],[rhs,rhs],'r-','LineWidth',1.5);
hold on;
for ii=1:100
    x0 = randn(n,1);

    % Fixed-Point
    [xfinal_FP_M1, x_iter_FP_M1, err_iter_FP_M1, runtime_FP_M1] = fixed_point(@q, M, x0,maxiter,tol);
    r_err_FP{ii} = max_n_to_end(r_linear_err(x_iter_FP_M1));
    
    % AA(m) 
    m_values = [1,2,3];
    for j=1:length(m_values)
        m = m_values(j);
        [xfinal_M1, x_iter_M1, err_iter_M1, runtime_M1] = AA_Rn(@q, M, x0, m, beta, maxiter, tol);
        r_err = max_n_to_end(r_linear_err(x_iter_M1));
        x_AA_iter{j} = x_iter_M1;
        r_err_AA{j,ii}  = r_err;
    end
    
    % Checking inequality in Eq. (4) 
    
    for k=1:length(m_values)
        x_iter_M1 = x_AA_iter{k};
        for i=2:length(x_iter_M1)-1
            num = norm((W-eye(n,n))*x_iter_M1{i+1});
            den = norm((W-eye(n,n))*x_iter_M1{i-1});
            lhs(i-1) = num/den;
        end
        LHS{k}=lhs;
        if rhs - max(lhs) < min_diff
            min_diff = rhs-max(lhs);
        end
        lhs = [];
    end
    plot(LHS{1},'b-','LineWidth',1.5);
    hold on;
    plot(LHS{2},'g-','LineWidth',1.5);
    hold on;
    plot(LHS{3},'m-','LineWidth',1.5);
    hold on;
end
legend({'RHS of Eq. (4)','LHS-AA(1)','LHS-AA(2)','LHS-AA(3)'},'FontSize',14,'FontWeight','bold','Location','best');
xlabel('Iterations (k)','FontSize',14,'FontWeight','bold');
ylabel('LHS(k) of Eq. (4)','FontSize',14,'FontWeight','bold')
axis([0,65 ,0, 0.70])
title('Tightness of Theorem 1','FontSize',16,'FontWeight','bold');

% Plot Results: R-Linear Convergence Rate
figure;
for ii=1:100
    plot(r_err_FP{ii},'b-','LineWidth',1.5);
    hold on;
    plot(r_err_AA{1,ii},'r-','LineWidth',1.5);
    hold on;
    plot(r_err_AA{2,ii},'k-','LineWidth',1.5);
    hold on;
    plot(r_err_AA{3,ii},'g-','LineWidth',1.5);
    hold on;
    plot([1,length(r_err_AA{1,11})],[r_est,r_est],'k--','LineWidth',1.5);
end
set(gca, 'YDir', 'reverse')
legend({'FP','AA(1)','AA(2)','AA(3)','Thm.1 Est.'},'FontSize',14,'FontWeight','bold','Location','best');
xlabel('Iterations (k)','FontSize',16,'FontWeight','bold');
ylabel('r_{est}(k)','FontSize',16,'FontWeight','bold')
axis([0,65 ,0.4, 1.8])
title('Comparing R-Linear Convergence Rates','FontSize',16,'FontWeight','bold');


%--------------------------------------------------------------------------
% Helper Fuctions
%--------------------------------------------------------------------------
function out = q(x,M)
    out = M*x;
end

function out = r_bound(a,maxW,minW)
    num = abs(maxW-minW)*a;
    den = sqrt(4+((2-maxW-minW)^2)*(a^2));
    term1 = asin(num/den);
    term2 = abs( atan(a) - atan((1-0.5*(maxW + minW))*a)  );
    out = term1 + term2;  
end

function [xfinal_FP, x_iter_FP, err_iter_FP, runtime_FP] = fixed_point(q,q_data,x0,max_iter,tol)
    tic;
    x = x0;
    xtemp = q(x0,q_data);
    err_iter_FP(1) = norm(x-xtemp);
    x = xtemp;
    x_iter_FP{1}=x0;
    x_iter_FP{2}=x;
    k = 1;
    while (k<=max_iter)&&(err_iter_FP(end)>=tol)
        xtemp = q(x,q_data);
        err_iter_FP(end+1) = norm(xtemp-x);
        x = xtemp;
        x_iter_FP{end+1}=x;
        k = k +1;
    end
    xfinal_FP = x;
    runtime_FP = toc;
    fprintf('\n_________________________________________\n')
    fprintf('_________________________________________\n')
    fprintf('Fixed-Point Iteration: \n')
    fprintf('\t Final Error Value: %7.4e\n',err_iter_FP(end))
    fprintf('\t Total Iterations: %5.0f\n',k-1)
    fprintf('\t Total time: %5.2f seconds \n',runtime_FP);
end

function r_err = r_linear_err(x_iter)
    for i=1:length(x_iter)
        r_err(i) = norm(x_iter{i})^(1/i);
    end
end

function new_vec = max_n_to_end(vec)
    new_vec = zeros(length(vec),1);
    for i=1:length(vec)
        new_vec(i) = max(vec(i:end));
    end
end
