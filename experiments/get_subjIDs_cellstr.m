function subjIDs = get_subjIDs_cellstr(Z)
nsamples = size(Z,1);
subjIDs = cell(nsamples,1);
for i =1:nsamples
    subjIDs{i} = num2str(find(Z(i,:) ==1));
end
