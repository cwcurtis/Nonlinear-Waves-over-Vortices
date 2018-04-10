function S = make_folder(rows,cols,mu,gam,F,tf,e)

kk = 1;

S = ['rows_' num2str(rows) '_cols_' num2str(cols) '_mu_pt' num2str(mu) '_gam_pt' num2str(gam) '_F_' num2str(F) '_tf_' num2str(tf) '_elp' num2str(e)];
S(regexp(S,'[.]'))=[];
temp = S;

while isequal(exist(temp,'dir'),7)
    temp = [S ' _(' num2str(kk) ')'];
    kk = kk+1;
end

S = temp;
mkdir(S)