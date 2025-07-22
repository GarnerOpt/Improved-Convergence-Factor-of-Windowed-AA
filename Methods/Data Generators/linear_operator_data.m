% This function forms the data for the linear function
% M1 - symmetric
% M2 - non-symmetric
function [M1,M2] = linear_operator_data(n,type,emax,emin)

    U2 = randn(n,n);
    U2 = U2*U2';

    A = randn(n,n);
    A = A*A';
    [U1,~]=eig(A);

    k=1;
    while k <= n
        temp = rand(1)-0.1;
        if abs(temp) < 0.9
            d(k,1)=temp;
            k = k+1;
        end
    end

    if type==1
        [~,min_i]=min(d);
        [~,max_i]=max(d);
        d(min_i,1)=emin;
        d(max_i,1)=emax;
    end
     d(end) = -0.95;    
     D = diag(d);
     M1 = U1*D*U1';
     M2 = U2*D*inv(U2);

    if type == 1
        M1 = D;
    end

end