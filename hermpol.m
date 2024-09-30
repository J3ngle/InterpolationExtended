function [herm,hxv]=hermpol(x,y,dy,z)
n=max(size(x));
m=max(size(z));
herm=[];
    if length(x) ~= length(y)
        error('Input vectors must have the same length.');
    end
for j=1:m
    xx=z(j);
    hxv=0;
        for i=1:n
            den=1;
            num=1;
            xn=x(i);
            derLi=0;
    for k=1:n
    if k~=i
        num=num*(xx-x(k));
        arg=xn-x(k);
        den=den*arg;
        derLi=derLi+1/arg;
    end
end
Lix2=(num/den).^2;
p=(1-2*(xx-xn)*derLi)*Lix2;
q=(xx-xn)*Lix2;
hxv=hxv+(y(i)*p+dy(i)*q);
end
herm=[herm,hxv];
end
return
