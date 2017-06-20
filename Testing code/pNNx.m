function pNNx1 = pNNx(RRI,x)
RR_diff = abs(diff(RRI));
num_diff_x = sum(RR_diff > x);
num_diff = numel(RR_diff);
pNNx1 = (num_diff_x/num_diff)*100;
end
