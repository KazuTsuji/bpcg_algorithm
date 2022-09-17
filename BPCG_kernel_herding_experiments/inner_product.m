%compute the inner product in an RKHS

function value=inner_product(c,d,p_c,p_d,epsilon,dim,functype)
if functype==1
    func = @(x,y,epsilon) exp(-epsilon*abs(x-y));
    d_function = @(x,y,epsilon) exp(-epsilon*norm(x-y,1));
elseif functype==2
    func = @(x,y,epsilon) exp(-epsilon*(x-y)^2);
    d_function = @(x,y,epsilon) exp(-epsilon*(norm(x-y))^2);
elseif functype==3
    d_function = @(x,y,epsilon) (1+norm(x-y))*exp(-norm(x-y));
elseif functype==4
    d_function = @(x,y,epsilon) (1+norm(x-y)+(norm(x-y) ^2 )/3)*exp(-norm(x-y));
elseif functype==5
    d_function = @(x,y,epsilon) exp(-epsilon*norm(x-y));
else
    'error_arises'
end
    m=length(c);
    n=length(d);
    value=0;
    for ii=1:m
        for jj=1:n
            value_ij=d_function(p_c(ii,:),p_d(jj,:),epsilon);
            value=value+c(ii)*d(jj)*value_ij;
        end
    end
end

