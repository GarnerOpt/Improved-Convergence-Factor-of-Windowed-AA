%==========================================================================
% Linear Operator - Verying Eq. (4) Tightness Test 2 
%==========================================================================
% GOAL: 
% Empirically verying the tightness of Eq.(4) for the case M is 5x5
% non-diagonal, symmetric. 
%--------------------------------------------------------------------------

fprintf('\n=============================================================\n')
fprintf('=============================================================\n')
fprintf('               LINEAR OPERATOR EXPERIMENT 1                 \n')
fprintf('=============================================================\n')
fprintf('=============================================================\n\n')

% Algorithm Settings 
maxiter = 1000; 
tol = 1e-12; 
m = 2;  
beta = 1;

% Data and randomized settings 
rng(13)
n = 5;
M = diag([-0.07   0.62    -0.55    -0.6   0.15  ]);  

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
sol_length = 0;

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
        sol_length = max([sol_length,length(LHS{k})]);
    end
    plot(LHS{1},'b-','LineWidth',1.5);
    hold on;
    plot(LHS{2},'g-','LineWidth',1.5);
    hold on;
    plot(LHS{3},'m-','LineWidth',1.5);
    hold on;
end
legend({'RHS of Eq. (4)','LHS-AA(1)','LHS-AA(2)','LHS-AA(3)'},'FontSize',14,'FontWeight','bold','Location','best');
xlabel('Iterations (k)','FontSize',16,'FontWeight','bold');
ylabel('Values in Eq. (4)','FontSize',16,'FontWeight','bold')
axis([0 sol_length+10 0 rhs+0.02])
title('Tightness of Theorem 1','FontSize',16,'FontWeight','bold');

fprintf('\n=============================================================\n')
fprintf('Minimum Difference Between LHS and RHS of Eq.(4): %5.4f\n',min_diff)
fprintf('Relative Percent of Minimum Gap: %3.2f\n',100*(min_diff/rhs))
fprintf('=============================================================\n')

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
