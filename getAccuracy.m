function [ out ] = getAccuracy(d)
% Calculate precision, recall, false alarm rate and f-measure between the first class vs the rest from multi class pr dataset
    [C, ~] = confmat(d);
    acc = sum(diag(C))/sum(sum(C));
    out = acc; 
end
