function explained_var = sseExplainedCal(signals,predicted_signals)
% signal = channel * tsample
% predicted_signal = channel*tsample
sse_residual = sum((signals - predicted_signals).^2,2);
sse_total = sum((signals-nanmean(signals,2)).^2,2);
explained_var = 1 - (sse_residual./sse_total);
end