q = 10; % number of observations
targets =  zeros(2,q);
for i = 1:q
    % if twin
   targets(1,i) = 1;
   % else
   targets(2,i) = 1;
end

outputs = rand(2,q);   % probability/thresholds, etc

[tpr,fpr,thresholds] = roc(targets,outputs)
plotroc(targets,outputs)