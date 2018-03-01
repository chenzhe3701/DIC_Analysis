
% summarize if tNote was used
% chenzhe, 2018-02-26

enable_used = 0;
disable_used = 0;
enabled = [0 0];
disabled = [0 0];
for ii=1:1417
   for jj = 1:length(stru(ii).c2t)
       if stru(ii).cEnable(jj) == 1
           enable_used = enable_used + 1;
           enabled = [enabled; stru(ii).gID, stru(ii).cLabel(jj)];
       end
       if stru(ii).cEnable(jj) == -1
           disable_used = disable_used + 1;
           disabled = [disabled; stru(ii).gID, stru(ii).cLabel(jj)];
       end
   end
end
display(enable_used);
display(disable_used);