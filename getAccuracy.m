function [ out ] = getAccuracy(d)
% Accuracy 
    [C, ~] = confmat(d);
    acc = sum(diag(C))/sum(sum(C));
    out = acc; 
end
