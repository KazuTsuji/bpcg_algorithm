%compute the worst case error

function error= derive_error(m,c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)

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
        
        error=0;
        %error=error+0.25*integral2(@(x,y) exp(-epsilon*(x-y).^2),-1,1,-1,1);
        if functype==1
            for jj=1:length(partion)
                for kk=1:length(partion)
                    val=1;
                    for dd=1:dim
                        val=val * 1/(c_var(jj,dd))* 1/(c_var(kk,dd))* integral2(@(x,y) (exp(-epsilon*abs(x-y))).*(exp(-var(jj)* (x-mean(jj,dd)).^2)).*exp(-var(kk)* (y-mean(kk,dd)).^2),-1,1,-1,1);
                    end
                    error=error+partion(jj)*partion(kk)*val;
                end
            end
        elseif functype==2
            for jj=1:length(partion)
                for kk=1:length(partion)
                    val=1;
                    for dd=1:dim
                        val=val * 1/(c_var(jj,dd))* 1/(c_var(kk,dd))* integral2(@(x,y) (exp(-epsilon*(x-y).^2)).*(exp(-var(jj)* (x-mean(jj,dd)).^2)).*exp(-var(kk)* (y-mean(kk,dd)).^2),-1,1,-1,1);
                    end
                    error=error+partion(jj)*partion(kk)*val;
                end
            end
        elseif functype==3
            if dim==2
                %error=load('matern32_2d.mat').value;
                error=load('matern32_2d_val.mat').matern32_2dval;
            elseif dim==3
                error=load('matern32_3d_val.mat').matern32_3dval;
                %error=load('matern32_3d.mat').value;
            end
        elseif functype==4
            if dim==2
                error=load('matern52_2d_integvalue.mat').value2d;
                %error=load('matern52_2d.mat').value;
            elseif dim==3
                error=load('matern52_3d.mat').value;
            end
        elseif functype==5
            if dim==2
                error=load('exp_2d.mat').value;
            elseif dim==3
                error=load('exp_3d.mat').value;
            end
            
        else
            'error_arises'
        end
        %ここはmuのノルム
        for i=1:m
            error=error-2*c(i)*mu(ind_num(i));
            for j=1:m
                value_ij=d_function(p(i,:),p(j,:),epsilon);
                error=error+c(i)*c(j)*value_ij;
            end
        end
        error=sqrt(error);
    end  
