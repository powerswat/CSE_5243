function [dist_mat] = calcDist(X, Y, order)

U=~isnan(Y); Y(~U)=0;
V=~isnan(X); X(~V)=0;
if order == 2
    dist_mat=sqrt(abs(X.^order*U' + V*Y'.^order - 2*X*Y'));
else
    dist_mat=abs(X*U' - V*Y');
end

end