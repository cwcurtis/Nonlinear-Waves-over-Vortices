function S = make_folder(rows,cols,mu,gam,F,tf)

kk = 1;

S = ['rows=' num2str(rows) ', cols=' num2str(cols) ', mu=' num2str(mu) ', gam=' num2str(gam) ', F=' num2str(F) ',tf=' num2str(tf)];

temp = S;
while isequal(exist(temp,'dir'),7)
    temp = [S ' (' num2str(kk) ')'];
    kk = kk+1;
end

S = temp;
mkdir(S)