b=unique(ind_mat(A_raw>0))
a=unique(ind_M)
LL = [];
for ii=1:92
    if ~ismember(b(ii),a)
        LL = [LL,b(ii)]
    end
end

%%

tp_ind(1,:)=ind_cell{1};
tp_ind(2,:)=ind_cell{2};

tp_skl(1,:)=skl_cell{1};
tp_skl(2,:)=skl_cell{2};