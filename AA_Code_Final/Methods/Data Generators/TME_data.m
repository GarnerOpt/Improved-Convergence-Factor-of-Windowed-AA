% This function forms the data for the TME numerical experiments 
% 
% Inputs (variable name - variable type - variable description)
% p - positive integer - dimension of the data points i.e. R^p
% n - positive integer - number of data points in R^p
% type - [1,2,3]       - selects the data generating model 


function [X,p,n] = TME_data(p,n,type)

    if type == 1
        X = zeros(p,n);
            for i=1:p
                for j=1:p
                    Sp(i,j) = (0.7)^(abs(i-j));
                 end
            end
            [U,D] = eig(Sp);
            Sp0_5 = U*(D.^(1/2))*U';

            for i=1:n
                X(:,i) = (Sp0_5)*randn(p,1);
            end

    % Inlier-Outlier Model 1        
    elseif type == 2
        n1=497;n0=500;D=100;d=50;
        p=D; n=n0+n1; %p=10, n=1995 
        X=[randn(n0,D)*randn(D,D);[rand(n1,d)/d^0.5,zeros(n1,D-d)]];
        X=X';


    % Inlier-Outlier Model 2 
    elseif type == 3
        p = 50;     % Dimension of the data i.e. \mathbb{R}^{400}
        n = 1120;   % Number of data points
        k=120;p0=5;
        X=[randn(k,p0),zeros(k,p-p0)]';
        X=[X,randn(p,n-k)];

    else
        fprintf('ERROR: This value of TYPE not supported. Choose TYPE=1,2 or 3. \n')
    end

end