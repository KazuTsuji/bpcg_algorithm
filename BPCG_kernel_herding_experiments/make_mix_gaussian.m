function make_mix_gaussian(dim,num,meanarr,partion,var)
    %{
    meanarr=randn(num,dim);
    partion=rand(1,num);
    var=30*rand(1,num);
    
    partion=partion/sum(partion)
    %}
    %{
    X=linspace(-1,1,100);
    Y=linspace(-1,1,100);
    Z=zeros(100,100);
    for ii=1:100
        for jj=1:100
            val=0;
            for kk=1:num
                val=val+partion(kk)*exp(-var(kk)*norm([X(ii),Y(jj)]-meanarr(kk,:))^2);
            end
            Z(ii,jj)=val;
        end
    end
    contourf(X,Y,Z,10,show)
    savefig('figures/mixture_gaussian.fig')
%}  
    function val=f(x,y)
        val=0;
        for kk=1:num
            val=val+partion(kk)*exp(-var(kk)*norm([x,y]-meanarr(kk,:))^2);
        end
    end
    fcontour(f,'Fill','on');
end