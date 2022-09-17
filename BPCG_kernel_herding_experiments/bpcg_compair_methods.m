%compair methods on the convergence speed of the worst case errors for the number of nodes when dim=2

%d=2
%
epsilon=1;
n=80;
dim=2;

meanarr=[[0,0]];
partion=[1];
var=[1];

%{
meanarr=[[0.8,0.3];[-0.2,-0.2];[-0.9,-0.9];[0.6,0.9]];
partion=[0.2,0.3,0.3,0.2];
var=[4,8,6,4];
%}

%make_mix_gaussian(dim,length(partion),meanarr,partion,var)
delta=0.0001;
K=20;
T=35;
functype=3;

option=4;
%}
maxtime=10000;

%
%d=3
%{
epsilon=1;
n=20;
dim=3;
meanarr=[[0,0,0]];
 
iteration=1;
 
partion=[1];
 
var=[1];
delta=0.0001;
K=20;
T=18;
functype=3;
option=4;

%}
%d=4
%{
epsilon=1;
n=10;
dim=4;
meanarr=[[0,0,0,0]];
 
iteration=1;
 
partion=[1];
 
var=[1];
delta=0.0001;
K=20;
T=18;
functype=2;
option=1;

%}
points =linspace(-1,1,n) ;
 
c_var=zeros(length(partion),dim);
for jj=1:dim
    for ii=1:length(partion)
        c_var(ii,jj)=1/sqrt(var(ii)) *(sqrt(pi)/2)*(erf(sqrt(var(ii))*(1-meanarr(ii,jj)))-  erf(sqrt(var(ii))*(-1-meanarr(ii,jj)))) ;   
    end
end
mu=zeros(1,n^(dim));
 
for ii=0:n^(dim)-1;
    indexes= func_index(ii,n,dim);
    value=1;    
    if functype==1;
        val=zeros(1,length(partion))+1;
        for jj=1:length(partion)
            for kk=1:dim
                int_value1= erf(sqrt(var(jj))* ( points(indexes(kk)) - (2*meanarr(jj,kk)*var(jj)+epsilon)/(2*var(jj) ) ) ) - erf(sqrt(var(jj))* ( -1 - (2*meanarr(jj,kk)*var(jj)+epsilon)/(2*var(jj) ) ) );
                int1= 1/(c_var(jj,kk)) * exp(-epsilon*points(indexes(kk)) -var(jj)* (meanarr(jj,kk))^2 )* exp((2*meanarr(jj,kk)*var(jj)+epsilon)^2 /(4*var(jj) )) * sqrt(pi)/2 *1/sqrt(var(jj))* int_value1;
                
                int_value2= erf(sqrt(var(jj))* ( 1 - (2*meanarr(jj,kk)*var(jj)- epsilon)/(2*var(jj) ) ) ) - erf(sqrt(var(jj))* ( points(indexes(kk)) - (2*meanarr(jj,kk)*var(jj)-epsilon)/(2*var(jj)) ) );
                int2= 1/(c_var(jj,kk)) * exp(epsilon*points(indexes(kk)) -var(jj)* (meanarr(jj,kk))^2 )* exp((2*meanarr(jj,kk)*var(jj)-epsilon)^2 /(4*var(jj) )) * sqrt(pi)/2 *1/sqrt(var(jj))* int_value2;
                
                val(jj)=val(jj)*(int1+int2);
            end
        end
        value= partion*transpose(val);
       
    elseif functype==2;
        val=zeros(1,length(partion))+1;
        for jj=1:length(partion)
            for kk=1:dim
                int_value= erf(sqrt(var(jj)+epsilon)+ (meanarr(jj,kk)*var(jj) +epsilon*points(indexes(kk)))/sqrt(var(jj)+epsilon) )+ erf(sqrt(var(jj)+epsilon) - (meanarr(jj,kk)*var(jj)+epsilon*points(indexes(kk)))/sqrt(var(jj)+epsilon) );
                val(jj)=val(jj)* 1/(c_var(jj,kk))* exp(-var(jj)*(meanarr(jj,kk))^2 - epsilon*(points(indexes(kk)))^2 + ((meanarr(jj,kk)*var(jj)+epsilon*(points(indexes(kk))))^2)/(var(jj)+epsilon) )* 1/sqrt(var(jj)+epsilon) * (sqrt(pi)/2)*int_value ; 
            end
        end
        value= partion*transpose(val);
    end
    mu(ii+1)=value;
    
    if functype==3;
        p_ar=[];
        for l=1:dim;
            p_ar=[p_ar,points(indexes(l))];
        end
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 ) ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu(ii+1)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 ) ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu(ii+1)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
    
    elseif functype==4;
        p_ar=[];
        for l=1:dim;
            p_ar=[p_ar,points(indexes(l))];
        end
        if dim==2
            ff=@(x,y)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2 )/3 ) .*exp(-sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu(ii+1)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z)(1+sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )+((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2 )/3 ) .*exp(-sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu(ii+1)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
        
    elseif functype==5;
        p_ar=[];
        for l=1:dim;
            p_ar=[p_ar,points(indexes(l))];
        end
        if dim==2
            ff=@(x,y) exp(-epsilon*sqrt( (x-p_ar(1)).^2+(y-p_ar(2)).^2 ));
            mu(ii+1)=(0.5)^dim *integral2(ff,-1,1,-1,1);
        elseif dim==3
            ff=@(x,y,z) exp(-epsilon*sqrt((x-p_ar(1)).^2+(y-p_ar(2)).^2+(z-p_ar(3)).^2));
            mu(ii+1)=(0.5)^dim *integral3(ff,-1,1,-1,1,-1,1);
        end
    end
    
end

order_1=[];
order_2=[];
order_45=[];
order_47=[];
num_nodes=[];


output_monte_carlo= monte_carlo(100,epsilon,dim,functype,option);
len_monte_carlo=  int64(length(output_monte_carlo)/2);
monte_carlo_err =[];
monte_carlo_nodes=[];

for i=1: len_monte_carlo
    m= int64(output_monte_carlo(i));
    monte_carlo_err = [monte_carlo_err, output_monte_carlo(len_monte_carlo +i)];
    monte_carlo_nodes=[monte_carlo_nodes,m];
    order_1=[order_1,1/double(m)];
    order_2=[order_2,1/sqrt(double(m))];
end


output_linesearch= linesearch(maxtime,100,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
len_ls=  int64(length(output_linesearch)/2);
line_search =[];
ls_nodes=[];

for i=1: len_ls
    m= output_linesearch(i);
    line_search = [line_search, output_linesearch(len_ls +i)];
    ls_nodes=[ls_nodes,m];
end

output_Away= Away(maxtime,100,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
len_away=  int64(length(output_Away)/2);
away =[];
away_nodes=[];

for i=1: len_away
    m= output_Away(i);
    away = [away, output_Away(len_away +i)];
    away_nodes=[away_nodes,m];
end

output_PWF= PWF(maxtime,100,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
len_pwf=  int64(length(output_PWF)/2);
pwf =[];
pwf_nodes=[];

for i=1: len_pwf
    m= output_PWF(i);
    pwf = [pwf, output_PWF(len_pwf +i)];
    pwf_nodes=[pwf_nodes,m];
end

output_eqweight= eqweight_herding(maxtime,100,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
len_eq=  int64(length(output_eqweight)/2);
eqweight =[];
eq_nodes=[];

for i=1: len_eq
    m= output_eqweight(i);
    eqweight = [eqweight, output_eqweight(len_eq +i)];
    eq_nodes=[eq_nodes,m];
end

output_SBQ= SBQ(maxtime,30,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
len_SBQ=  int64(length(output_SBQ)/2);
SBQ_err =[];
SBQ_nodes=[];
order_exp=[];

for i=2: len_SBQ
    m= output_SBQ(i);
    SBQ_err = [SBQ_err, output_SBQ(len_SBQ +i)];
    SBQ_nodes=[SBQ_nodes,m];
end
%}

order_45=[];
order_47=[];
order_exp=[];
order_nodes=[];

for jj=1:35
    m=double(jj);
    order_nodes=[order_nodes,m];
    order_45=[order_45,(1/double(m))^(5/4)];
    order_47=[order_47,(1/double(m))^(7/4)];
    order_exp=[order_exp,exp(-double(m)^(1/dim))];
end


output_bcg_pairwise_linesearch= bcg_pairwise_linesearch(maxtime,2000,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
len_bcg_pairwise_ls=  int64(length(output_bcg_pairwise_linesearch)/2);
bcg_pairwise_line_search =[];
bcg_pairwise_ls_nodes=[];

for i=1: len_bcg_pairwise_ls
    m= output_bcg_pairwise_linesearch(i);
    bcg_pairwise_line_search = [bcg_pairwise_line_search, output_bcg_pairwise_linesearch(len_bcg_pairwise_ls +i)];
    bcg_pairwise_ls_nodes=[bcg_pairwise_ls_nodes,m];
end

output_bcg_pairwise_linesearch_lazified= bcg_pairwise_linesearch_lazified(maxtime,5000,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
len_bcg_pairwise_ls_lazified=  int64(length(output_bcg_pairwise_linesearch_lazified)/2);
bcg_pairwise_line_search_lazified =[];
bcg_pairwise_ls_nodes_lazified=[];

for i=1: len_bcg_pairwise_ls_lazified
    m= output_bcg_pairwise_linesearch_lazified(i);
    bcg_pairwise_line_search_lazified = [bcg_pairwise_line_search_lazified, output_bcg_pairwise_linesearch_lazified(len_bcg_pairwise_ls_lazified +i)];
    bcg_pairwise_ls_nodes_lazified=[bcg_pairwise_ls_nodes_lazified,m];
end
 

%newcolors = {'#F00','#ff8c00','#FF0','#0B0','#00F','#00ffff','#000','#0072BD','#A2142F'};
%colororder(newcolors)

semilogy(monte_carlo_nodes,monte_carlo_err,'Color','#7E2F8E','LineWidth',2);
hold on

semilogy(eq_nodes,eqweight,'-','Color','#00F','LineWidth',2);
hold on

semilogy(ls_nodes,line_search,'-','Color','#d2691e','LineWidth',2);
hold on

semilogy(away_nodes,away,'-','Color','#ff69b4','LineWidth',2);
hold on

semilogy(pwf_nodes,pwf,'-','Color','#0B0','LineWidth',2);
hold on

semilogy(SBQ_nodes,SBQ_err,'-','Color','#00ffff','LineWidth',2);
hold on

semilogy(bcg_pairwise_ls_nodes,bcg_pairwise_line_search,'Color','#000','LineWidth',2);
hold on
semilogy(bcg_pairwise_ls_nodes_lazified,bcg_pairwise_line_search_lazified,'Color','#F00','LineWidth',2);
hold on

%semilogy(monte_carlo_nodes,order_2,'Color','#A2142F')
%hold on

semilogy(order_nodes,order_45,'--','Color','#A2142F')
hold on
%semilogy(order_nodes,order_47,'--','Color','#A2142F')
%hold on
%semilogy(order_nodes,order_exp,'--','Color','#A2142F')
%hold on

xlabel('the number of nodes','FontSize',20)
%xlabel('the number of iterations','FontSize',20)
%xlabel('computational time (s)','FontSize',20)
ylabel('MMD','FontSize',20)
 
hold off

legend({'Monte carlo','eq-weight','linesearch','Away','PCG','SBQ','BPCG','BPCG-lazified','(1/n)^{5/4}'},'FontSize',18,'NumColumns',2)
%legend({'Monte carlo','eq-weight','linesearch','Away','PCG','SBQ','BPCG','BPCG-lazified','(1/n)^{7/4}'},'FontSize',14,'NumColumns',2)
%legend({'eq-weight','linesearch','Away','PCG','SBQ','BPCG','BPCG-lazified','exp(-n^{1/2})'},'FontSize',18,'NumColumns',2)
%legend({'eq-weight','linesearch','Away','PCG','SBQ','BPCG','BPCG-lazified'},'FontSize',20,'NumColumns',2)
%legend({'eq-weight','linesearch','Away','PCG','BPCG','BPCG-lazified'},'FontSize',20,'NumColumns',2)