
% chenzhe, 2018-01-31
% Based on matlab's algorithm using kmeans++ to select initial guess
% points, with part of them initially provided.

function c_init = kmeans_pp_init(X,k,nPages,c_provided)

if ~isreal(X)
    error(message('stats:kmeans:ComplexData'));
end
wasnan = any(isnan(X),2);
hadNaNs = any(wasnan);
if hadNaNs
%     warning(message('stats:kmeans:MissingDataRemoved'));
    X = X(~wasnan,:);
end

if(isempty(X))
   warning('empty sample space, randomize.');
   X = rand(k,size(c_provided,2));
end

if size(c_provided,1)>1
c_provided = c_provided(1:k,:);
warning('too many initial guess points provided. select only k of them.');
end

[n,p] = size(X);

index = zeros(1,k);
c0(1:size(c_provided,1),:) = c_provided;
c0(size(c_provided,1)+1:k,:) = 0;
minDist = inf(n,1);

for ip = 1:nPages
    % Select the rest of the seeds by a probabilistic model
    for ii = size(c_provided,1)+1:k
        minDist = min(minDist,distfun(X,c0(ii-1,:),'sqeuclidean'));
        denominator = sum(minDist);
        if denominator==0 || isinf(denominator) || isnan(denominator)
            c0(ii:k,:) = datasample(X,k-ii+1,1,'Replace',false);
            break;
        end
        sampleProbability = minDist/denominator;
        [c0(ii,:), index(ii)] = datasample(X,1,1,'Replace',false,'Weights',sampleProbability);
    end
    c_init(:,:,ip) = c0;
end


end

%------------------------------------------------------------------

function D = distfun(X, C, dist, iter,rep, reps)
%DISTFUN Calculate point to cluster centroid distances.

switch dist
    case 'sqeuclidean'
        D = pdist2(X,C,'seuclidean');  
%     case {'cosine','correlation'}
%         % The points are normalized, centroids are not, so normalize them
%         normC = sqrt(sum(C.^2, 1));
%         if any(normC < eps(class(normC))) % small relative to unit-length data points
%             if reps==1
%                 error(message('stats:kmeans:ZeroCentroid', iter));
%             else
%                 error(message('stats:kmeans:ZeroCentroidRep', iter, rep));
%             end
%             
%         end
%         C = bsxfun(@rdivide,C,normC);
%         D = pdist2mex(X,C,'cos',[],[],[]);  
end
end % function