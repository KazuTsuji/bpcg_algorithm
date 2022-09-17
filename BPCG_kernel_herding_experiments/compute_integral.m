function value= compute_integral(func,dim,m)
    %等間隔点で数値積分構成
    value=0;
    for ii=1:m
        point1=rand(1,dim)*2 -1;
        point2=rand(1,dim)*2 -1;
        value=value+func(point1,point2);
    end
    value=value/m;
end        