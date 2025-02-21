function Score = IdentityMatchFactor(s1,s2)
C=999;

m2 = sum(s1~=0 & s2~=0); %number of non zero elements in s1 and s2

matchFactor = C*dot(s1,s2)/(sqrt(sum(s1.^2,"all"))*sqrt(sum(s2.^2,"all")));
%calculate composite score
ds1=[0,s1(1:end-1)];
ds2=[0,s2(1:end-1)];
R = (s1.*ds1).*(s2.*ds2);
gamma1=s1./ds1.*ds2./s2;
gamma1(~isfinite(gamma1))=0; %remove NaN and Inf values
R(gamma1==0)=0;
R(R>0) = gamma1(R>0);
alpha = R~=0;
m1 = sum(alpha); %number of nonzero elements in R
R = R(alpha);
if isempty(R)
    F = 0;
else
    F = (sum(R*min([R,1./R],[],"all"),'all'))/sum(R,"all");
end
Score = C*(m1*F+m2*matchFactor/C)/(m1+m2);
Score = round(Score);
end