function S = make_folder(K,mu,gam,tf)

kk = 1;

S = ['_K_' num2str(K) '_mu_pt' num2str(mu) '_gam_pt' num2str(gam) '_tf_' num2str(tf)];
S(regexp(S,'[.]'))=[];
temp = S;

while isequal(exist(temp,'dir'),7)
    temp = [S '_(' num2str(kk) ')'];
    kk = kk+1;
end

S = temp;
mkdir(S)