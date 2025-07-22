%==========================================================================
% Anderson Acceleration Algorithm for fixed-point problems w/ restarting
%==========================================================================

%--------------------------------------------------------------------------
% Description:
%--------------------------------------------------------------------------
% This algorithm implements Anderson Acceleration with restarts. The algorithm is 
% set-up to apply to fixed point operators q mapping vectors from R^n to R^n. 

%--------------------------------------------------------------------------
% Inputs: (format: input name - input type - input descriptor)
%--------------------------------------------------------------------------
%  x0      - vector in R^n                - int. vector for algorithm 
%  q       - function mapping R^n to R^n  - operator/fixed point map
%  q_data  - cell_array                   - information to define q function
%  m       - positive integer             - depth of the algorithm until restarting 
%  beta    - [0,1]                        - dampening parameter 
%  maxiter - positive integer             - max. number of algorithm
%  tol     - positive real number         - termination tolerence of alg.

%--------------------------------------------------------------------------
% Outputs:
%--------------------------------------------------------------------------
% xfinal   - vector in R^n   - final iterate of alg. 
% err_iter - vector          - error of each iter. of alg, i.e. ||x-q(x)||
% runtime  - real number     - run time (in secs.) of alg. 

%==========================================================================

function [xfinal, x_iter, err_iter, runtime] = AArestart_Rn(q, q_data, x0, m, beta, maxiter, tol)
    tic; 
    k = 1;
    err_iter(1)=inf;
    x1   = q(x0,q_data);
    X{k} = x0;
    while (k <= maxiter) && (err_iter(end) >= tol)
        %--------------------------------------------------
        % Step 1. - First iteration via fixed-point method
        %--------------------------------------------------
        X{end+1} = x1;
        qX{1} = x1;
        if k==1
            err_iter(1) = norm(qX{1}-X{1},2);
            x_iter{1}=x0;
            x_iter{2}=x1;
        else
            err_iter(end+1)= norm(qX{1}-X{1},2);
            x_iter{end+1}=x1;
        end
    
        %-------------------------------
        % Step 2. - Main Loop of Alg. 1
        %-------------------------------
        k_in = 1;
        while (k <= maxiter) && (err_iter(end) >= tol) && (k_in <= m)
            qX{end+1} = q(X{end},q_data);
            err_iter(end+1) = norm(qX{end}-X{end},2);    
            Residuals = form_residuals(X, qX);
            [R,r] = construct_residual_matrix(Residuals);
            alpha = anderson_step(R,r);
            newPoint  = new_anderson_point(X,qX,alpha,beta);
            X{end+1}  = newPoint;
            x_iter{end+1}=newPoint;

             if length(qX) == m+1
                qX = {};
                X  = {};
            end

            k = k+1;
            k_in = k_in + 1;
        end
        x0 = x_iter{end};
        x1 = q(x0,q_data);
        X{1} = x0;
    end
        %----------------------------------------------------
        % Step 3. - Store and print final iterate information
        %----------------------------------------------------
        xfinal = newPoint;
        runtime = toc;
        fprintf('_________________________________________\n')
        fprintf('Anderson Acceleration (m=%2.0f):\n',m)
        fprintf('\t Final Error Value: %7.4e\n',err_iter(end))
        fprintf('\t Total Iterations: %5.0f\n',length(x_iter))
        fprintf('\t Total time: %5.2f seconds \n',runtime)
end

%--------------------------------------------------------------------------
% Helper Functions
%--------------------------------------------------------------------------

function [alpha] = anderson_step(R,r) 
    [~,k]=size(R);
    alpha_short  = (R'*R+0*eye(k,k))\(-R'*r);
    alpha(1:k)   = alpha_short;
    alpha(end+1) = 1 - ones(k,1)'*alpha_short;
end

function newPoint = new_anderson_point(X,qX,alpha,beta)
    [p,~]=size(X{1});
    newPoint = zeros(p,1);
    for i=1:length(X)
        newPoint = newPoint + (1-beta)*alpha(i)*X{i} + beta*alpha(i)*qX{i};
    end
end

function Residuals = form_residuals(X, qX)
    for i=1:length(qX)
        Residuals{i} = qX{i} - X{i};
    end
end

function [R,r] = construct_residual_matrix(Residuals)
    n = length(Residuals{1});
    m = length(Residuals);
    R = zeros(n,m-1);
    r = Residuals{end};
    for j = 1:m-1
        R(:,j) = Residuals{j}-r;
    end
end

function new_array = slice_array(Xarray)
    new_array = cell(1,length(Xarray)-1);
    for i=2:length(Xarray)
        new_array{i-1} = Xarray{i};
    end
end


