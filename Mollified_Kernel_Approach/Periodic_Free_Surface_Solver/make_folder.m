function S = make_folder(rows,cols,K,mu,gam,F,tf,modu,kap,zoffc)

kk = 1;

S = ['rows_' num2str(rows) '_cols_' num2str(cols) '_K_' num2str(K) '_mu_pt' num2str(mu) '_gam_pt' num2str(gam) '_F_' num2str(F) '_tf_' num2str(tf) '_modu_' num2str(modu) '_kap_' num2str(kap) '_zoffc_' num2str(zoffc)];
S(regexp(S,'[.]'))=[];
temp = S;

while isequal(exist(temp,'dir'),7)
    temp = [S '_(' num2str(kk) ')'];
    kk = kk+1;
end

S = temp;
mkdir(S)