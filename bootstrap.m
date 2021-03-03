function [CI_L, CI_U, M] = bootstrap( data, rep, bins)
    fprintf('bootstrapping ... \n'); tic
    x = size(data,2); y = size(data,1); m = zeros(rep, x); x2 = 0:x-1;
    if(isempty(bins)==1); bins = ones(1, x)*y; end
    
    for b=1:rep % choose random sets of data
        fprintf('\b\b\b\b\b%3.0f%%\n', b/rep*100)
        m(b,:) = mean(data(ceil(rand(y,x).*bins) + x2*y),1);
    end
    m(isnan(m)==1)=0;
    for xi = 1:x % sort data through time points
        m(:,xi) = sort(m(:,xi));
    end
    m = m(ceil(rep*0.05)+1:floor(rep-rep*0.05),:); % cut of top and bottom 5%
    CI_L = m(1,:); CI_U = m(end,:); M = mean(m,1);
    fprintf('\b\b\b\b\bcompleted in %.0f seconds\n', toc)
    
end

