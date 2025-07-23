nback_idxs = [];
for row = 1:length(data.trialinfo)

    if row > 1

        if ((data.trialinfo(row,2) == data.trialinfo(row-1,2)) && (data.trialinfo(row,1) == data.trialinfo(row-1,1)))

            if data.trialinfo(row,3) == data.trialinfo(row-1,3)
                nback_idxs = [nback_idxs row];
            end
        end
    end
end
