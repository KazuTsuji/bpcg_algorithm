function output = SBQ(Maxtime,m,epsilon,n,dim,points,mu,c_var,mean,var,partion,functype,option)
tstart=tic;
t_calc=0.0;
time=[];

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

error_ar=[];
number_of_points=[];

syms x;


p=zeros(m,dim);%array of nodes
ind=zeros(m,dim);%array of indexes of nodes
ind_num=[];
c=[];%array of weights

same_count=0;

nodes=[];
error_nodes=[];
iteration_t=[];
error_t=[];
mu_norm=(derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype))^2
i=1;

time_total=0;
while ( (length(c)+same_count) < m)&&(time_total <Maxtime );
    if i==1;
        a=find(mu==max(mu));
        ind_num=[ind_num,a(1)];
        number=a(1)-1;
       
        for jj=1:dim;
            ind(i,jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
            p(i,jj)=points(ind(i,jj));
            number= rem(number,n^(dim-jj));
        end
  
    else;
        z_t_array=zeros(1,length(mu));
        z_matrix=zeros(length(mu),length(ind_num)+1);
        for jj=0:n^(dim)-1;
            if  (numel(find(ind_num == jj+1))==0 );
                z=zeros(1,length(ind_num)+1);
                K_matrix=zeros(length(ind_num)+1,length(ind_num)+1);

                x_point=zeros(1,dim);
                numb=jj;
                for tt=1:dim
                    x_ind=idivide(int64(numb),int64(n^(dim-tt)) )+1;
                    x_point(tt)=points(x_ind);
                    numb= rem(numb,n^(dim-tt));
                end

                K_matrix(length(ind_num)+1,length(ind_num)+1)=d_function(x_point,x_point,epsilon);
                z(length(ind_num)+1)=mu(jj+1);
                for kk=1:length(ind_num)
                    K_matrix(length(ind_num)+1,kk)=d_function(x_point,p(kk,:),epsilon);
                    K_matrix(kk,length(ind_num)+1)=d_function(x_point,p(kk,:),epsilon);
                    z(kk)=mu(ind_num(kk));
                    for ll=1:length(ind_num)
                        K_matrix(kk,ll)=d_function(p(ll,:),p(kk,:),epsilon);
                    end
                end
                if cond(K_matrix) < 10^10 
                    z_t_array(jj+1)= mu_norm - z*inv(K_matrix)*transpose(z);
                    z_matrix(jj+1,:)=transpose(inv(K_matrix)*transpose(z));
                else
                    z_t_array(jj+1)=10000;
                end
            else
                z_t_array(jj+1)=10000;
            end
        end

        a=find(z_t_array==min(z_t_array));
        min(z_t_array);
        ind_num_i = a(1);
        p_i=zeros(1,dim);
        ind_i=zeros(1,dim);
        number=a(1)-1;
        
        for jj=1:dim;
            ind_i(jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
            p_i(jj)=points(ind_i(jj));
            number= rem(number,n^(dim-jj));
        end
    end
    
    if i==1;
        alpha=mu(ind_num(1));
        c=[c,alpha];
        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    else
        c=z_matrix(ind_num_i,:);
        ind_num=[ind_num,ind_num_i];
        p(length(c),:)=p_i;
        ind(length(c),:)=ind_i;

        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    end
    
    iteration_t=[iteration_t,i];
    error_t=[error_t, derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    
    time=[time,double((toc(tstart)-t_calc))];
    tic
    t_calc=double(t_calc)+double(toc);
    
    error_ar=[error_ar,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    number_of_points=[number_of_points,length(c)+sum(same_count)];
    
    
    i=i+1;
    time_total=double((toc(tstart)-t_calc));
end

if option==1
    output=[nodes,error_nodes];
elseif option==2
    output=[iteration_t,error_t];
elseif option==3
    output=[time,error_t];
elseif option==4
    output=[number_of_points, error_ar];
end

end