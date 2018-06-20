function embedding = diffusionMap(mDist,args)

sigma               = args.eps * median(mDist(:));
aff                 = exp(-mDist.^2/(2*sigma^2)); 
aff                 = stochastic(aff);
eigsnum             = args.eigsnum;
[eigvecs, eigvals]  = eigs(aff, eigsnum);
embedding           = eigvecs(:,2:eigsnum)*eigvals(2:eigsnum,2:eigsnum);
