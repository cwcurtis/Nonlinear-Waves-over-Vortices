function S = make_folder(rows,cols,K,mu,gam,F,tf,Gamma,coords,pattern)

kk = 1;

if strcmp(coords, 'polar')
    S = ['polar__rows=' num2str(rows) ', cols=' num2str(cols) ', K=' num2str(K) ', mu=' num2str(mu) ', gam=' num2str(gam) ', F=' num2str(F) ',tf=' num2str(tf) ',GAMMA=' num2str(Gamma) ',pattern=' num2str(pattern)];
else
    S = ['rows=' num2str(rows) ', cols=' num2str(cols) ', K=' num2str(K) ', mu=' num2str(mu) ', gam=' num2str(gam) ', F=' num2str(F) ',tf=' num2str(tf) ',GAMMA=' num2str(Gamma) ',pattern=' num2str(pattern)];
end

temp = S;
while isequal(exist(temp,'dir'),7)
    temp = [S ' (' num2str(kk) ')'];
    kk = kk+1;
end

S = temp;
mkdir(S)