function [ aff_mat ] = CalcInitAff( data, params )
% Calculates the affinity between the columns of the data as thresholded 
% cosine similarity between columns, or as a Gaussian kernel of width 
% eps*(median distance between 5 nearest neighbors of all points).
% 
% Inputs:
%     data - M-by-N matrix whose columns are N points in R^M 
%     params - struct array with user parameters  
% 
% Output:
%     aff_mat - N-by-N symmetric matrix of non-negative affinities 
%--------------------------------------------------------------------------

if strcmp(params.metric,'cosine_similarity'),
    inner_products = data.'*data;
    [ij_norm, ji_norm] = meshgrid(diag(inner_products));
    aff_mat = inner_products./sqrt(ij_norm.*ji_norm);
    aff_mat(aff_mat < params.thresh) = 0;
else
    data = data.';
    euc_dist = squareform(pdist(data));
    [~, nn_dist] = knnsearch(data, data, 'k', params.knn);
    sigma = params.eps * median(nn_dist(:));
    aff_mat = exp(-euc_dist.^2/(2*sigma^2));    
end

% figure, imagesc(aff_mat), colormap jet, colorbar, axis on, hold on
% if params.on_rows,
%     title('Initial Row Affinity'), hold off
% else
%     title('Initial Col Affinity'), hold off
% end

end

