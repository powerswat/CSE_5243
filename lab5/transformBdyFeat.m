function [bdyFeatMat, bdyVectLabel] = transformBdyFeat(bdyVectLabel, ...
                                        bdy_is, bdy_js, bdy_vals)

% Convert the topic and body lists to matrices
bdyFeatMat = zeros(max(bdy_is),length(bdyVectLabel));
for i=1:length(bdy_js)
    if exist('bdy_vals','var')
        bdyFeatMat(bdy_is(i),bdy_js(i)) = bdy_vals(i);
    else
        bdyFeatMat(bdy_is(i),bdy_js(i)) = 1;
    end
end

end