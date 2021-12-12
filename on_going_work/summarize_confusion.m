function tbl = summarize_confusion(confMat)
% input confusion matrix = [tp,fp,fn,tn], do summarize
% chenzhe, 2019-12-26

t = confMat;
[nR,nC] = size(t);
t(nR+1,:) = sum(t,1);
acc = (t(:,1)+t(:,4))./sum(t,2);
ppr = t(:,1)./(t(:,1)+t(:,2));  % precision
tpr = t(:,1)./(t(:,1)+t(:,3));  % sensitivity

tbl = array2table([t,acc,ppr,tpr]);
tbl.Properties.VariableNames = {'TP','FP','FN','TN','Accuracy','Precision','Sensitivity'};
tbl.Properties.RowNames{nR+1} = 'Total';

disp(tbl);

end