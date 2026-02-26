function [vxRaw,vyRaw] = HS_flowfield(tracePhase,useGPU)
%% frameN in the first dimension
frameN = size(tracePhase,1);

vxRaw = zeros(frameN-1,size(tracePhase,2),size(tracePhase,3));
vyRaw = zeros(frameN-1,size(tracePhase,2),size(tracePhase,3));

if useGPU
    tracePhase = gpuArray(tracePhase);
    vxRaw = gpuArray(vxRaw);
    vyRaw = gpuArray(vyRaw);
end

for k = 1:(frameN-1)
    tic
    A1 = squeeze(tracePhase(k,:,:));
    A2 = squeeze(tracePhase(k+1,:,:));
    [vxRaw(k,:,:), vyRaw(k,:,:)] = HS_phase_mod(A1, A2);
    if mod(k,100) == 0
        fprintf('frame %g/%g; time elapsed %g seconds \n', [k, frameN, toc])
    end
end
if useGPU
    vxRaw = gather(vxRaw);
    vyRaw = gather(vyRaw);
    tracePhase = gather(tracePhase);
end