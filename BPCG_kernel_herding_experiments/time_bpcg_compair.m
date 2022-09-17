%compute the convergence for computational time
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

option=3;
%}
maxtime=60;
max_time_copy=maxtime;
maxnodes=100000;
max_nodes_copy=maxnodes;
time_num=10;

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


stride=10;
for tt=1:time_num
    output_linesearch= linesearch(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
    len_ls=  int64(length(output_linesearch)/2);
    line_search =[];
    ls_nodes=[];

    fact=0;
    for i=1: len_ls
        if fact*stride < i
            line_search = [line_search, output_linesearch(len_ls +i)];
            ls_nodes=[ls_nodes,output_linesearch(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        line_search_ave=line_search;
        ls_nodes_ave=ls_nodes;
        maxtime=10000;
        maxnodes=len_ls;
    else
        line_search_ave=line_search_ave+line_search;
        ls_nodes_ave=ls_nodes_ave+ls_nodes;
    end
end
line_search=line_search_ave/time_num;
ls_nodes=ls_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;

for tt=1:time_num
    output_eqweight= eqweight_herding(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
    len_eq=  int64(length(output_eqweight)/2);
    eqweight =[];
    eq_nodes=[];

    fact=0;
    for i=1: len_eq
        if fact*stride < i
            eqweight = [eqweight, output_eqweight(len_eq +i)];
            eq_nodes=[eq_nodes,output_eqweight(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        eqweight_ave=eqweight;
        eq_nodes_ave=eq_nodes;
        maxtime=10000;
        maxnodes=len_eq;
    else
        eqweight_ave=eqweight_ave+eqweight;
        eq_nodes_ave=eq_nodes_ave+eq_nodes;
    end
end
eqweight=eqweight_ave/time_num;
eq_nodes=eq_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;

for tt=1:time_num
    output_Away= Away(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
    len_away=  int64(length(output_Away)/2);
    away =[];
    away_nodes=[];

    fact=0;
    for i=1: len_away
        if fact*stride < i
            away = [away, output_Away(len_away +i)];
            away_nodes=[away_nodes,output_Away(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        away_ave=away;
        away_nodes_ave=away_nodes;
        maxtime=10000;
        maxnodes=len_away;
    else
        away_ave=away_ave+away;
        away_nodes_ave=away_nodes_ave+away_nodes;
    end
end
away=away_ave/time_num;
away_nodes=away_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;



stride=10;

for tt=1:time_num
    output_pfw= PWF(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
    len_pfw=  int64(length(output_pfw)/2);
    pfw =[];
    pfw_nodes=[];

    fact=0;
    for i=1: len_pfw
        if fact*stride < i
            pfw = [pfw, output_pfw(len_pfw +i)];
            pfw_nodes=[pfw_nodes,output_pfw(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        pfw_ave=pfw;
        pfw_nodes_ave=pfw_nodes;
        maxtime=10000;
        maxnodes=len_pfw;
    else
        pfw_ave=pfw_ave+pfw;
        pfw_nodes_ave=pfw_nodes_ave+pfw_nodes;
    end
end
pfw=pfw_ave/time_num;
pfw_nodes=pfw_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;



stride=10;

for tt=1:time_num
    output_bcg_pw_lazified= bcg_pairwise_linesearch_lazified(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
    len_bcg_pw_lazified=  int64(length(output_bcg_pw_lazified)/2);
    bcg_pw_lazified =[];
    bcg_pw_lazified_nodes=[];

    fact=0;
    for i=1: len_bcg_pw_lazified
        if fact*stride < i
            bcg_pw_lazified = [bcg_pw_lazified, output_bcg_pw_lazified(len_bcg_pw_lazified +i)];
            bcg_pw_lazified_nodes=[bcg_pw_lazified_nodes,output_bcg_pw_lazified(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        bcg_pw_lazified_ave=bcg_pw_lazified;
        bcg_pw_lazified_nodes_ave=bcg_pw_lazified_nodes;
        maxtime=10000;
        maxnodes=len_bcg_pw_lazified;
    else
        bcg_pw_lazified_ave=bcg_pw_lazified_ave+bcg_pw_lazified;
        bcg_pw_lazified_nodes_ave=bcg_pw_lazified_nodes_ave+bcg_pw_lazified_nodes;
    end
end
bcg_pw_lazified=bcg_pw_lazified_ave/time_num;
bcg_pw_lazified_nodes=bcg_pw_lazified_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;


stride=10;

for tt=1:time_num
    output_bcg_pw= bcg_pairwise_linesearch(maxtime,maxnodes,epsilon,n,dim,points,c_var,mu,var,meanarr,partion,functype,option);
    len_bcg_pw=  int64(length(output_bcg_pw)/2);
    bcg_pw =[];
    bcg_pw_nodes=[];

    fact=0;
    for i=1: len_bcg_pw
        if fact*stride < i
            bcg_pw = [bcg_pw, output_bcg_pw(len_bcg_pw +i)];
            bcg_pw_nodes=[bcg_pw_nodes,output_bcg_pw(i)];
            fact=fact+1;
        end
    end
    
    if tt==1
        bcg_pw_ave=bcg_pw;
        bcg_pw_nodes_ave=bcg_pw_nodes;
        maxtime=10000;
        maxnodes=len_bcg_pw;
    else
        bcg_pw_ave=bcg_pw_ave+bcg_pw;
        bcg_pw_nodes_ave=bcg_pw_nodes_ave+bcg_pw_nodes;
    end
end
bcg_pw=bcg_pw_ave/time_num;
bcg_pw_nodes=bcg_pw_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;

stride=1;
for tt=1:time_num
    output_SBQ=  SBQ(maxtime,maxnodes,epsilon,n,dim,points,mu,c_var,meanarr,var,partion,functype,option);
    len_SBQ=  int64(length(output_SBQ)/2);
    SBQ_err =[];
    SBQ_nodes=[];

    fact=0;
    for i=1: len_SBQ
        if fact*stride < i
            SBQ_err = [SBQ_err, output_SBQ(len_SBQ +i)];
            SBQ_nodes=[SBQ_nodes,output_SBQ(i)];
            fact=fact+1;
        end
    end
    if tt==1
        SBQ_ave=SBQ_err;
        SBQ_nodes_ave=SBQ_nodes;
        maxtime=10000;
        maxnodes=len_SBQ;
    else
        SBQ_ave=SBQ_ave+SBQ_err;
        SBQ_nodes_ave=SBQ_nodes_ave+SBQ_nodes;
    end
end
SBQ_err=SBQ_ave/time_num;
SBQ_nodes=SBQ_nodes_ave/time_num;
maxtime=max_time_copy;
maxnodes=max_nodes_copy;

semilogy(eq_nodes,eqweight,'Color','#00F','LineWidth',2);
hold on

semilogy(ls_nodes,line_search,'-','Color','#d2691e','LineWidth',2);
hold on

semilogy(away_nodes,away,'-','Color','#ff69b4','LineWidth',2);
hold on

semilogy(pfw_nodes,pfw,'-','Color','#0B0','LineWidth',2);
hold on

semilogy(SBQ_nodes,SBQ_err,'-','Color','#00ffff','LineWidth',2);
hold on

semilogy(bcg_pw_nodes,bcg_pw,'Color','#000','LineWidth',2);
hold on
semilogy(bcg_pw_lazified_nodes,bcg_pw_lazified,'Color','#F00','LineWidth',2);
hold on


xlabel('computational time (s)','FontSize',20)
ylabel('MMD','FontSize',20)

legend({'eq-weight','linesearch','Away','PCG','SBQ','BPCG','BPCG-lazified'},'FontSize',20,'NumColumns',2)
