function rank_es = AIC(eigvalue,bet)

L = 1/(1-bet);
AIC_temp = 1e6;

n = length(eigvalue);
rank_es = 1;
for k = 0 : n-1
    
    sum_eig_noise = sum(eigvalue(k+1:end));
    ari_mean = sum_eig_noise/(n-k);
    
    prod_eig_noise = prod(eigvalue(k+1:end));
    geo_mean = prod_eig_noise^(1/(n-k));
    ak = ari_mean/geo_mean;
    
    AIC_val = (n-k)*L*log(ak) + k*(2*n-k); 
    if AIC_val  <  AIC_temp
        AIC_temp = AIC_val;
        rank_es = k;
    end
end

end