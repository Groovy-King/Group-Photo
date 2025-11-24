function [ res ] = func( x )

for k=1 :length(x)
if x(k)>1 
        res(k)=1-((2/sqrt(x(k)^2-1))*atan(sqrt((x(k)-1)/(x(k)+1))));
else if x(k)<1
        res(k)=1-((2/sqrt(1-x(k)^2))*atanh(sqrt((1-x(k))/(1+x(k)))));
    else
        res(k)=0;
    end
end
end

end

