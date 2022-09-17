%kernel herding with linesearch

function output = bcg_pairwise_linesearch_lazified(Maxtime,T,epsilon,n,dim,points,c_var,mu,var,mean,partion,functype,option)

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

syms x;

number_of_points=[];


p=zeros(T,dim);%array of nodes
ind=zeros(T,dim);%array of indexes of nodes
ind_num=[];
c=[];%array of weights
i=1;
same_count=[];

drop_count=0;
%memorize the inner products and norms at the previous iteration of t
x_t_norm=0;
mu_x_t =0;
x_t_array=zeros(1,length(mu));
nn=1;
trunc=10^(-6);

nodes=[];
error_nodes=[];
iteration_t=[];
error_t=[];


K_accr=1;

total_time=0;
t=1;
while (t < (T+1))&&(total_time< Maxtime)
    if i==1;
        a=find(mu==max(mu));
        ind_num=[ind_num,a(1)];
        number=a(1)-1;
        
        for jj=1:dim;
            ind(i,jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
            p(i,jj)=points(ind(i,jj));
            number= rem(number,n^(dim-jj));
        end
        
    end
    %compute the step size
    v1=0;
    v2=0;
    if i==1;
        alpha=1;
        c=[c,alpha];
        same_count=[same_count,0];
        x_t_norm= 1; %norm |x_t|
        mu_x_t= mu(a(1));%<mu,x_t>
        add_number=length(c);
        i=i+1;
        nodes=[nodes,length(c)];
        error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
        
        mu_prime=mu;
        for j=0:n^(dim)-1;
            indexes=func_index(j,n,dim);
            points_ar=[];
            for l=1:dim;
                points_ar=[points_ar,points(indexes(l))];
            end
            value=d_function(points_ar,p(1,:),epsilon);
            mu_prime(j+1)=mu_prime(j+1) - value;
        end
        
        Phi= ( x_t_norm - mu_x_t + max(mu_prime) )/2;
    else
        x_t_norm=0;
        mu_x_t=0;
        
        x_t_array_St=zeros(1,length(c));
        for jj=1:length(c);
            mu_x_t=mu_x_t+c(jj)*mu(ind_num(jj));
            for kk=1:length(c)
                value=d_function(p(jj,:),p(kk,:),epsilon);
                x_t_norm=x_t_norm+c(kk)*c(jj)*value;
                x_t_array_St(jj)=x_t_array_St(jj)+ c(kk)*value;
            end
        end 
        
        b_St=x_t_array_St-mu(ind_num);
        
        st=find(b_St==min(b_St));%fw_direction
        at=find(b_St==max(b_St));%away_direction
        
        fw_value=   b_St(st(1));
        away_value=  b_St(at(1));
        
        if  ( (away_value - fw_value) > Phi) 
            %move in pairwise direction%
            fw_aw_value=d_function(p(at(1),:),p(st(1),:),epsilon);
            
            val_frac= 2-2*fw_aw_value;
            
            alpha= (away_value - fw_value) /val_frac;
            lambda= c(at(1));
       
            if (alpha < lambda ) 
                c(st(1))=c(st(1))+alpha;
                c(at(1))=c(at(1))-alpha;
            else
                c(st(1))=c(st(1))+lambda;
                c(at(1))=[];
                c= c/sum(c);
                same_count(at(1))=[];
                p(at(1),:)=[];
                ind_num(at(1))=[];
                ind(at(1),:)=[];
                drop_count= drop_count+1;
            end
            
        else
            x_t_array=zeros(1,length(mu));
            for j=0:n^(dim)-1;
                indexes=func_index(j,n,dim);
                for kk=1:length(c)
                    points_ar=[];
                    for l=1:dim;
                        points_ar=[points_ar,points(indexes(l))];
                    end
                    value=d_function(points_ar,p(kk,:),epsilon);
                    x_t_array(j+1)=x_t_array(j+1)+ c(kk)*value;
                end
            end
            mu_x_t=0;
            x_t_norm=0;
            
            
            for jj=1:length(c)
                mu_x_t=mu_x_t+c(jj)*mu(ind_num(jj));
                for kk=1:length(c)
                    value=d_function(p(jj,:),p(kk,:),epsilon);
                    x_t_norm=x_t_norm+c(jj)*c(kk)*value;
                end
            end
            
            b=x_t_array - mu;
            wt=find(b==min(b));
            ind_num_wt = wt(1);
            p_wt=zeros(1,dim);
            ind_wt=zeros(1,dim);
            number=wt(1)-1;
            
            for jj=1:dim;
                ind_wt(jj)=idivide(int64(number),int64(n^(dim-jj)) )+1;
                p_wt(jj)=points(ind_wt(jj));
                number= rem(number,n^(dim-jj));
            end
            
            mu_v = mu(ind_num_wt);%<mu,v>
            x_t_v= inner_product(c,[1],p,p_wt,epsilon,dim,functype); %<x_t,v>
        
        
            v1=x_t_norm -mu_x_t - x_t_v +mu_v ;%<x_t -mu, x_t -v_t>
            v2=x_t_norm - 2*x_t_v  +1;%|x_t -v_t|^2
            alpha=min([1,(v1/v2)]);
            
            
            if (v1 > Phi/K_accr)
                x_t_norm= (1-alpha)^2 *x_t_norm +2*(1-alpha)*alpha*x_t_v + alpha^2;%|(1-alpha)x_t +alpha v|^2
                mu_x_t= (1-alpha)*mu_x_t + alpha*mu_v; 
        
      
                bbb =find(ind_num == ind_num_wt);
                if (numel(bbb)~=0)
                    c= (1-alpha)*c;
                    c(bbb(1))=c(bbb(1))+alpha;
                    same_count(bbb(1))=same_count(bbb(1))+1;
                    add_number=bbb(1);
                else
                    c= [(1-alpha)*c,alpha];
                    ind_num= [ind_num, ind_num_wt];
                    same_count=[same_count,0];
                    for jj=1:dim
                        ind(length(c),jj)=ind_wt(jj);
                        p(length(c),jj)=points(ind_wt(jj));
                        add_number=length(c);
                    end
                end
            
                nodes=[nodes,length(c)];
                error_nodes=[error_nodes,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
            
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
            
            else
                Phi=Phi/2;
            end
        end
    end

    error_ar=[error_ar,derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    number_of_points=[number_of_points,length(c)+sum(same_count)];
    
    iteration_t=[iteration_t,t];
    error_t=[error_t, derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)];
    time=[time,double((toc(tstart)-t_calc))];
    tic
    t_calc=double(t_calc)+double(toc);
    
    total_time=double((toc(tstart)-t_calc));
    t=t+1;
end

derive_error(length(c),c,ind_num,dim,p,epsilon,c_var,mean,var,mu,partion,functype)
sum(same_count)
if option==1
    output=[nodes,error_nodes];
elseif option==2
    output=[iteration_t,error_t];
elseif option==3
    output=[time,error_t];
elseif option==4
    output=[number_of_points, error_ar];
end

drop_count

end