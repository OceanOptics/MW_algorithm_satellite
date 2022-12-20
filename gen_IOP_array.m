function [IOP_matrix,a_sw,dim] = gen_IOP_array(rrs,U,std,bandwidth,SNR,nm,temp,S,Y,ap443, ap750,bbp700)
% selecting spectral ranges 
rrs          = rrs(:,(nm>=630 & nm <=670 | nm>=700  & nm <=2500));
U            = U(:,(nm>=630 & nm <=670 | nm>=700  & nm <=2500));
std          = std(:,(nm>=630 & nm <=670 | nm>=700  & nm <=2500));
bandwidth    = bandwidth(nm>=630 & nm <=670 | nm>=700 & nm <=2500);
SNR          = SNR(nm>=630 & nm <=670 | nm>=700 & nm <=2500);
nm           = nm(nm>=630 & nm <=670 | nm>=700 & nm <=2500);

% number of combination of IOPs
dim = (length(S)*length(Y)*length(ap443)*length(ap750)*length(bbp700)); 

%generate vectors for IOP's -> water absorption:

[absorp_cor] = asw_correction(temp,nm,bandwidth);
a_sw = repmat(absorp_cor',1,1,dim);
a_sw = permute(a_sw, [1 3 2]);
clear absorp_cor

%-------------------------------------------------------------------------%
%generate vectors for IOP's -> absorption:
for i = 1:length(S)
    for ii = 1:length(nm)
        for iii = 1:length(ap443)
            for iiii = 1:length(ap750)
                a_p_sat(ii,iii,iiii,i)      =  ((ap443(iii) .* (exp( -S(i) * ((nm(ii)+bandwidth(ii)./2)-443)))) -  (ap443(iii) .* (exp(-S(i) * ((nm(ii)-bandwidth(ii)./2)-443)))) ) ./ (-S(i).*bandwidth(ii));
                a_p_offset(ii,iii,iiii,i)   =   (ap443(iii) .* (- (exp(-S(i)*(750-443)))) )  + ap750(iiii);
                S_slope(ii,iii,iiii,i) = S(i);
                ap443_matrix(ii,iii,iiii,i) = ap443(iii);
                ap750_matrix(ii,iii,iiii,i) = ap750(iiii);
            end
        end
    end
end

a_p = a_p_sat + a_p_offset;
 
ap = reshape(a_p,[length(nm),(length(S)*length(ap443)*length(ap750))])';
S_res = reshape(S_slope,[length(nm),(length(S)*length(ap443)*length(ap750))])';
a_nap443_res = reshape(ap443_matrix,[length(nm),(length(S)*length(ap443)*length(ap750))])';
a_nap750_res = reshape(ap750_matrix,[length(nm),(length(S)*length(ap443)*length(ap750))])';
clear i ii iii a_p

%-------------------------------------------------------------------------%
%generate vectors for IOP's -> backscattering:
for n = 1:length(Y)
    for nn = 1:length(nm)
        for nnn = 1:length(bbp700)
            if Y ~= 1
                bb_p(nn,nnn,n) = (bbp700(nnn).*(((nm(nn)+bandwidth(nn)./2)./700).^(-Y(n)+1))  -  bbp700(nnn).*(((nm(nn)-bandwidth(nn)./2)./700).^(-Y(n)+1)))./((-Y(n) + 1).*bandwidth(nn));
                Y_slope(nn,nnn,n) = Y(n);
            else
                bb_p(nn,nnn,n) = (bbp700(nnn).*700.* (log((nm(nn)+bandwidth(nn)./2) ./(nm(nn)-bandwidth(nn)./2)))) ./ bandwidth(nn);
                Y_slope(nn,nnn,n) = Y(n);
            end
            bbp_matrix(nn,nnn,n) = bbp700(nnn);
        end
    end
end

bbp = reshape(bb_p,[length(nm),(length(Y)*length(bbp700))])';
Y_res = reshape(Y_slope,[length(nm),(length(Y)*length(bbp700))])';
bbp_res = reshape(bbp_matrix,[length(nm),(length(Y)*length(bbp700))])';
clear bb_p n nn nnn


%-------------------------------------------------------------------------%
%generate all possible combination of eigenvectors
k = 0;
for i = 1:length(S)*length(ap443)*length(ap750)
    for n = 1:length(Y)*length(bbp700)
        k = k+1;
        IOP_matrix(:, :, k) = [ap(i, :); bbp(n, :)];
        shape_matrix(:, :, k) = [S_res(i, :); Y_res(n, :)];
        coefs_matrix(:, :, k) = [a_nap443_res(i, :); a_nap750_res(i, :); bbp_res(n, :)];
    end
end

end

