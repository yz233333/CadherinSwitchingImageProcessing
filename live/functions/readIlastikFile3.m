function [mask_1,mask_2,mask_3] = readIlastikFile3(filename,numLabel3)
% read ilastik h5 file and output masks of label 1, 2(, and 3)
% Input: filename; numLabel3 - if 3 labels, numLabel3=1, if 2 labels, =0
% PS: label 1=membrane, label 2=background, label 3=non-membrane
% cell
% PS, select dead cells: label 1=non-background, label 2=background 

immask = h5read(filename, '/exported_data');
immask = squeeze(immask);

nd = ndims(immask);
perm_array = 1:nd;
perm_array(1:2) = [2 1];

mask_1 = immask == 1;
mask_1 = permute(mask_1,perm_array);
mask_2 = immask == 2;
mask_2 = permute(mask_2,perm_array);
if numLabel3
    mask_3 = immask == 3;
    mask_3 = permute(mask_3,perm_array);
end
end