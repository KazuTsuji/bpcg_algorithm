%kernel herding with linesearch

function output = Away(T,m,epsilon,n,dim,points,c_var,mu,var,mean,partion,functype,option)
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
FW_AW_flag="FW";

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
        same_count=[same_count,0];
        c=[c,alpha];
        x_t_norm= 1; %norm |x_t|
        mu_x_t= mu(a(1));%<mu,x_t>
        add_number=length(c);
        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
        FW_AW_flag="FW";
    else
        v_fw= -min(b) -mu_x_t +x_t_norm;
        b_St=b(ind_num);
        v_aw= max(b_St) -x_t_norm+mu_x_t;
        if v_aw > v_fw
            FW_AW_flag="AWAY";
            x_t_array_St=x_t_array(ind_num);
            ind_num_aw = find(b_St==max(b_St));
            v1=v_aw;
            v2=x_t_norm-2*x_t_array_St(ind_num_aw)+d_function(p(ind_num_aw,:),p(ind_num_aw,:),epsilon);
            alpha=min([ c(ind_num_aw)/(1-c(ind_num_aw)),(v1/v2)]);
            x_t_norm= (1+alpha)^2 *x_t_norm -2*(1+alpha)*alpha*x_t_array_St(ind_num_aw)+ alpha^2 *d_function(p(ind_num_aw,:),p(ind_num_aw,:),epsilon);%|(1+alpha)x_t -alpha v|^2
            mu_St=mu(ind_num);
            mu_x_t= (1+alpha)*mu_x_t - alpha*mu_St(ind_num_aw);
            if alpha < (c(ind_num_aw)/(1-c(ind_num_aw)))
                c=(1+alpha)*c;
                c(ind_num_aw)=c(ind_num_aw)-alpha;
            else
                c=(1+alpha)*c;
                c(ind_num_aw)=[];
                c= c/sum(c);
                same_count(ind_num_aw)=[];
                p(ind_num_aw,:)=[];
                ind_num(ind_num_aw)=[];
                ind(ind_num_aw,:)=[];
            end
        else
            FW_AW_flag="FW";
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

            mu_v = mu(ind_num_i);%<mu,v>
            x_t_v= inner_product(c,[1],p,p_i,epsilon,dim,functype); %<x_t,v>
        
        
            v1=x_t_norm -mu_x_t - x_t_v +mu_v ;%<x_t -mu, x_t -v_t>
            v2=x_t_norm - 2*x_t_v  +1;%|x_t -v_t|^2
            alpha=min([1,(v1/v2)]);
        
            x_t_norm= (1-alpha)^2 *x_t_norm +2*(1-alpha)*alpha*x_t_v + alpha^2;%|(1-alpha)x_t +alpha v|^2
            mu_x_t= (1-alpha)*mu_x_t + alpha*mu_v; 
        
            bbb =find(ind_num == ind_num_i);
            if(numel(bbb)~=0)
                c= (1-alpha)*c;
                c(bbb(1))=c(bbb(1))+alpha;
                same_count(bbb(1))=same_count(bbb(1))+1;
                add_number=bbb(1);
            else
                c= [(1-alpha)*c,alpha];
                same_count=[same_count,0];
                ind_num= [ind_num, ind_num_i];
                for jj=1:dim
                    ind(length(c),jj)=ind_i(jj);
                    p(length(c),jj)=points(ind_i(jj));
                    add_number=length(c);
                end
            end
        end
        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,sqrt(mu_norm -2*mu_x_t + x_t_norm)];
    end
    if FW_AW_flag=="FW"
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
    elseif FW_AW_flag=="AWAY"
        x_t_array=(1+alpha)*x_t_array;
        for j=0:n^(dim)-1;
            indexes=func_index(j,n,dim);
            points_ar=[];
            for l=1:dim;
                points_ar=[points_ar,points(indexes(l))];
            end
            value=d_function(points_ar,p(ind_num_aw,:),epsilon);
            x_t_array(j+1)=x_t_array(j+1)- alpha*value;
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