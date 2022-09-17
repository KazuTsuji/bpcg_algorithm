%kernel herding with equal weights

function output = monte_carlo(m,epsilon,dim,functype,option)


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
c=[];%array of weights

nodes=[];
error_nodes=[];
iteration_t=[];
error_t=[];

density_const= (1/ (erf(1)*sqrt(pi)) )^dim;

mu_ar=zeros(1,m);

if functype==1;
   val=1;
   for dd=1:dim
        val=val * integral2(@(x,y) exp(-epsilon*abs(x-y)).* exp(- x.^2).*exp(-y.^2),-1,1,-1,1);
    end
    mu_norm=val*density_const^2;
elseif functype==2;
    val=1;
    for dd=1:dim
        val=val *integral2(@(x,y) exp(-epsilon*(x-y).^2).*exp(- x.^2).*exp(- y.^2),-1,1,-1,1);
    end
    mu_norm=val*density_const^2
elseif functype==3
    if dim==2
        mu_norm=load('matern32_2d_val.mat').matern32_2dval;
    elseif dim==3
        mu_norm=load('matern32_3d_val.mat').matern32_3dval;
    end
elseif functype==4
    if dim==2
        mu_norm=load('matern52_2d_integvalue.mat').value2d;
    elseif dim==3
        mu_norm=load('matern52_3d.mat').value;
    end
elseif functype==5
    if dim==2
        mu_norm=load('exp_2d.mat').value;
    elseif dim==3
        mu_norm=load('exp_3d.mat').value;
    end
end

tstart=tic;
t_calc=0.0;
time=[];

mu_val_total=0;
quad_total=0;

for ii=1:m
    if functype==1 || functype==2
        flag=false;
        while flag==false
            point=randn(1,dim);
            ran_flag=true;
            for jj=1:dim
                if (point(jj) >1 )|| (point(jj)<-1 )
                    ran_flag=false;
                end
            end
            if ran_flag==true
                flag=true;
            end
        end
    elseif functype==3 || functype==4 || functype==5
        point=2*rand(1,dim)-1;
    end

    p(ii,:)=point;

    if functype==1;
        val=1;
        for kk=1:dim
            val=val*int(exp(-x^2)*exp(-epsilon*abs(x-point(kk))), -1, 1) ;
        end
        mu_ar(ii)=val*density_const;
       
    elseif functype==2;
        val=1;
        for kk=1:dim
            val=val*int(exp(-x^2)*exp(-epsilon*(x-point(kk))^2), -1, 1) ;
        end
        mu_ar(ii)=val*density_const;
    
    elseif functype==3;
        p_ar=point;
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 ) ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu_ar(ii)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 ) ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu_ar(ii)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
    
    elseif functype==4;
        p_ar=point;
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2 )/3 ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu_ar(ii)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )/3 ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu_ar(ii)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
        
    elseif functype==5;
        p_ar=point;
        if dim==2
            ff=@(x,y) exp(-epsilon*sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu_ar(ii)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z) exp(-epsilon*sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu_ar(ii)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
    end

    mu_val_total=mu_val_total+mu_ar(ii);
    quad_total=quad_total+d_function(p(ii,:),p(ii,:),epsilon);
    if ii>1
        for jj=1:ii-1
            quad_total= quad_total+ 2*d_function(p(jj,:),p(ii,:),epsilon);
        end
    end

    error_value=mu_norm - 2* 1/ii * mu_val_total + (1/ii)^2 *quad_total;
    error_value=sqrt(error_value);

    nodes=[nodes,ii];
    error_nodes=[error_nodes,error_value];
    
    iteration_t=[iteration_t,i];
    error_t=[error_t, error_value];
    
    time=[time,double((toc(tstart)-t_calc))];
    tic
    t_calc=double(t_calc)+double(toc);
    
    error_ar=[error_ar,error_value];
    number_of_points=[number_of_points,ii];

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

