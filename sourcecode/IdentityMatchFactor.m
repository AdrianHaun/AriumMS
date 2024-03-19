function Score = IdentityMatchFactor(s1,s2)
C=999;

m2 = sum(s1~=0 & s2~=0); %number of non zero elements in s1 and s2

matchFactor = C*dot(s1,s2)/(sqrt(sum(s1.^2,"all"))*sqrt(sum(s2.^2,"all")));
%calculate composite score
ds1=[0,s1(1:end-1)];
ds2=[0,s2(1:end-1)];
R1 = s1.*ds1.*s2.*ds2;
gamma1=s1./ds1.*ds2./s2;
gamma1(R1==0)=[];
if isempty(gamma1)
    gamma1 = 1;
end
m1 = numel(gamma1); %number of elements in gamma
F1 = (sum(gamma1*min([gamma1,1./gamma1],[],"all"),'all'))/sum(gamma1,"all");
Score = C*(m1*F1+m2*matchFactor/C)/(m1+m2);
end