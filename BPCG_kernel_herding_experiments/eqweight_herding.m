%kernel herding with equal weights

function output = eqweight_herding(T,m,epsilon,n,dim,points,mu,c_var,mean,var,partion,functype,option)
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


x_t_array=zeros(1,length(mu));
i=1;
time_total=0;
while ( (length(c)+same_count) < m)&&(time_total < T)
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
        x_t_array=(1-alpha)*x_t_array;
        for j=0:n^(dim)-1;
            indexes=func_index(j,n,dim);
            points_ar=[];
            for l=1:dim;
                points_ar=[points_ar,points(indexes(l))];
            end
            value=d_function(points_ar,p(add_number,:),epsilon);
            x_t_array(j+1)=x_t_array(j+1)+ alpha*value;
        end
        b=x_t_array- mu;
        
        a=find(b==min(b));
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
        alpha=1;
        c=[c,alpha];
        add_number=length(c);
    else
        kk=double(i);
        alpha=1/kk;
    
        bbb =find(ind_num == ind_num_i);
       
        if (numel(bbb)~=0)
            c= (1-alpha)*c;
            c(bbb(1))=c(bbb(1))+alpha;
            same_count=same_count+1;
            add_number=bbb(1);
        else
            c= [(1-alpha)*c,alpha];
            ind_num= [ind_num, ind_num_i];
            for jj=1:dim
                ind(length(c),jj)=ind_i(jj);
                p(length(c),jj)=points(ind_i(jj));
            end
            add_number=length(c);
            nodes=[nodes,length(c)];
            error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
        end
    
    end
    
    iteration_t=[iteration_t,i];
    error_t=[error_t, derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    
    time=[time,double((toc(tstart)-t_calc))];
    tic
    t_calc=double(t_calc)+double(toc);
    
    error_ar=[error_ar,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    number_of_points=[number_of_points,length(c)+sum(same_count)];
    
    
    time_total=double((toc(tstart)-t_calc));
    i=i+1;
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

