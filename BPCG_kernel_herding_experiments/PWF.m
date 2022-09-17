%kernel herding with linesearch

function output = PWF(T,m,epsilon,n,dim,points,c_var,mu,var,mean,partion,functype,option)
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

syms x;

error_ar=[];
number_of_points=[];


p=zeros(m,dim);%array of nodes
ind=zeros(m,dim);%array of indexes of nodes
ind_num=[];
c=[];%array of weights
i=1;
same_count=[];

mu_norm=(derive_error(0,[],[],dim,[],epsilon,c_var,mean,var,mu,partion,functype))^2;%|mu|_k ^2

%memorize the inner products and norms at the previous iteration of t
x_t_norm=0;
mu_x_t =0;
x_t_array=zeros(1,length(mu));

nodes=[];
error_nodes=[];
iteration_t=[];
error_t=[];

time_total=0;
while ( i < (m+1))&&(time_total < T)
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
        b=x_t_array- mu;
    end
    %compute the step size
    v1=0;
    v2=0;
    if i==1;
        alpha=1;
        c=[c,alpha];
        same_count=[0];
        x_t_norm= 1; %norm |x_t|
        mu_x_t= mu(a(1));%<mu,x_t>
        add_number=length(c);
        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    else
        %compute FW
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

        %compute FW
        x_t_array_St=x_t_array(ind_num);
        b_St=b(ind_num);
        aa = find(b_St==max(b_St));
        ind_num_aw=aa(1);
        v1= b_St(ind_num_aw)-b(ind_num_i);
        v2= d_function(p(ind_num_aw,:),p(ind_num_aw,:),epsilon)+d_function(p_i,p_i,epsilon)-2*d_function(p(ind_num_aw,:),p_i,epsilon);
        alpha=min([c(ind_num_aw),(v1/v2)]);
        
        x_t_norm= x_t_norm +2*alpha*(x_t_array(ind_num_i)-x_t_array_St(ind_num_aw))+alpha^2 *v2;
        mu_St=mu(ind_num);
        mu_x_t= mu_x_t + alpha*(mu(ind_num_i)-mu_St(ind_num_aw));

        %AWの更新
        if alpha < c(ind_num_aw)
            c(ind_num_aw)=c(ind_num_aw)-alpha;
            d_flag=false;
        else
            d_flag=true;
        end
        
        bbb =find(ind_num == ind_num_i);
        if (numel(bbb)~=0)
            ind_point=bbb(1);
            c(ind_point)=c(ind_point)+alpha;
            same_count(ind_point)=same_count(ind_point)+1;
            add_number=ind_point;
        else
            c= [c,alpha];
            ind_num= [ind_num, ind_num_i];
            same_count=[same_count,0];
            for jj=1:dim
                ind(length(c),jj)=ind_i(jj);
                p(length(c),jj)=points(ind_i(jj));
                add_number=length(c);
            end
        end
        if d_flag==true
            p_ind_num_aw=p(ind_num_aw,:);
            if add_number > ind_num_aw
                add_number=add_number-1;
            end
            c(ind_num_aw)=[];
            same_count(ind_num_aw)=[];
            p(ind_num_aw,:)=[];
            ind_num(ind_num_aw)=[];
            ind(ind_num_aw,:)=[];
        end

        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    end
    
    if i==1
        for j=0:n^(dim)-1;
            indexes=func_index(j,n,dim);
            points_ar=[];
            for l=1:dim;
                points_ar=[points_ar,points(indexes(l))];
            end
            value=d_function(points_ar,p(add_number,:),epsilon);
            x_t_array(j+1)=value;
        end
    else
        for j=0:n^(dim)-1;
            indexes=func_index(j,n,dim);
            points_ar=[];
            for l=1:dim;
                points_ar=[points_ar,points(indexes(l))];
            end
            value1=d_function(points_ar,p(add_number,:),epsilon);
            if d_flag==true
                value2=d_function(points_ar,p_ind_num_aw,epsilon);
            else
                value2=d_function(points_ar,p(ind_num_aw,:),epsilon);
            end
            x_t_array(j+1)=x_t_array(j+1)+alpha*(value1-value2) ;
        end
    end

    iteration_t=[iteration_t,i];
    error_t=[error_t, sqrt(mu_norm -2*mu_x_t + x_t_norm)];
    
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