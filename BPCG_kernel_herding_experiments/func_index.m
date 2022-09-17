%compute the indexes of a point
function indexes= func_index(number,n,dim)
    indexes=[];
    for jj=1:dim;
        indexes=[indexes,idivide(int64(number),int64(n^(dim-jj)) )+1];
        number= rem(number,n^(dim-jj));
    end