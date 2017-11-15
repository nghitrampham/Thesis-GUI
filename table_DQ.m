get_group3 = clusts_STD;

group1_idx=get_group3{1};

gender = {'female'; 'male'; 'female'; 'male'};

group  = [1;1;2;2];

female_group1_idx = group1_idx(group1_idx <= 63);
male_group1_idx = group1_idx(group1_idx >63);
DQ_female1_mean = nanmean(DQ(female_group1_idx));
DQ_male1_mean = nanmean(DQ(male_group1_idx));

group2_idx=get_group3{2};

female_group2_idx = group2_idx(group2_idx <= 63);
male_group2_idx = group2_idx(group2_idx >63);
DQ_female2_mean = nanmean(DQ(female_group2_idx));
DQ_male2_mean = nanmean(DQ(male_group2_idx));

 size_per_group = [length(female_group1_idx);length(male_group1_idx);length(female_group2_idx);...
            length(male_group2_idx)];       

DQ_list = [DQ_female1_mean;DQ_male1_mean;DQ_female2_mean;DQ_male2_mean];

group3_normal_GTable_sigma06 = table(group, gender, DQ_list, size_per_group);

% writetable(group3_normal_GTable_sigma06,'group3_normal_GTable_sigma06.csv')