A = dataset('XLSFile','status.csv');
C = dataset('XLSFile','merged_pvc.csv');
joinDiff = join(A,C,'LeftKeys',{'RID'},'RightKeys',{'RID'},'Type','fullouter','MergeKeys',true);
export(joinDiff,'XLSFile','merged.csv')
