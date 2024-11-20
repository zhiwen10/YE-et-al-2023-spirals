% get raw data from a set of U(usually a selected subspace from original U),V
% and compute a new set of Unew and Vnew
% useful to use the Vnew for reduced dimension regression analysis,etc.

function [Unew,Vnew,Sv,totalVar] = redoSVD(U,V)
% first average nearby nt0 (17 for 84592frames) frames to get downsized raw data and a new U
% because widefield data size is huge to fit into memory
% and do bacth processing
NavgFramesSVD = 5000;
ntotframes = size(V,2);
nt0 = ceil(ntotframes / NavgFramesSVD);
NavgFramesSVD = floor(ntotframes/nt0);
imgbatchSize = nt0 * floor(1000/nt0);
nbatch = floor(ntotframes/imgbatchSize);
%% read U, V and average nearby frames 
uSize = size(U,1);
ix = 0;
mov = zeros(uSize , NavgFramesSVD, 'single');
tic
for ibatch = 1:nbatch
        fprintf(1, '   frame %d out of %d\n', ix*nt0, ntotframes);
        batchStart = 1+(ibatch-1)*imgbatchSize;
        batchEnd = ibatch*imgbatchSize;
        data = U*V(:,batchStart:batchEnd);      
%         irange = 1:nt0*floor(size(data,2)/nt0);
%         data = data(:, irange);
        data = reshape(data, uSize, nt0, []);
        davg = single(squeeze(mean(data,2)));
        movRange = ix + (1:size(davg,2));
        mov(:,movRange) = davg;
        ix = ix + size(davg,2);
end
toc
mov(:, :, (ix+1):end) = [];
%% get COV matrix for mov
% nSVD = 200;
mov = mov-mean(mov,2);
COV             = mov' * mov/size(mov,1);
% total variance of data. If you ask for all Svs back then you will see
% this is equal to sum(Sv). In this case Sv are the singular values *of the
% covariance matrix* not of the original data - they are equal to the Sv of
% the original data squared (the variances per dimension). 
totalVar = sum(diag(COV));                             
% nSVD = min(size(COV,1)-2, nSVD);
%% svd for COV matrix to get new Ua
nSVD = 2000;
useGPU = 0;
if nSVD<1000 || size(COV,1)>1e4
    [Va, Sv]          = eigs(double(COV), nSVD);
else
    if useGPU
        [Va, Sv]         = svd(gpuArray(double(COV)));
        Va = gather(Va);
        Sv = gather(Sv);
    else
         [Va, Sv]         = svd(COV);
    end
    Va               = Va(:, 1:nSVD);
    Sv              = Sv(1:nSVD, 1:nSVD);
end
%
clear COV
Ua               = normc(mov * Va);
clear mov
Ua              = single(Ua);
Sv              = single(diag(Sv));
%% get new V (Fs) from new U and data, batch by batch
ix = 0;
Fs = zeros(nSVD, ntotframes, 'single');   

for ibatch = 1:nbatch+1    
    fprintf(1, '   frame %d out of %d\n', ix, ntotframes);  
    batchStart = 1+(ibatch-1)*imgbatchSize;
    if ibatch == nbatch+1 & mod(ntotframes,imgbatchSize)
        batchEnd = ntotframes;
    else
        batchEnd = ibatch*imgbatchSize;
    end
    data = U*V(:,batchStart:batchEnd);    

    % subtract mean as we did before
    data = data-mean(data,2);
    FsRange = ix + (1:size(data,2));
    Fs(:, FsRange) = Ua' * data;       
    ix = ix + size(data,2);
end
%%
Vnew = Fs; clear Fs
Unew = Ua;
end