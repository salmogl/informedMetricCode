function [ embedding ] = calcEmbedding( dist , args )

eucDist    = squareform(dist);
eps        = args.eps*median(dist);
aff        = exp(-eucDist.^2/eps);

D          = diag(sum(aff,2));
aff_stoch  = D\aff;

eigsnum            = 3;
[eigvecs, eigvals] = eigs(aff_stoch, eigsnum);
embedding          = eigvecs(:,2:3)*eigvals(2:3,2:3);

end

