clear
figure(1)
clf
%
% blasius equation
% f''' + 1/2 f f'' = 0
%
% b.c. f(0)=f'(0)=0 and f'(inf)=1
% 
% g = f'
% g' = f'' = h
% g'' = f''' = h' = -1/2 f f'' = -1/2 f h
%
etavec(1)=0;
%
deta=0.001;
total=10;
%
fvec(1)=0;
gvec(1)=0;
%
% vary initial value of h
%
hvec(1)=0.3319;
%
for i=1:total/deta;
%
etavec(i+1)=etavec(i)+deta;
%
fvec(i+1)=fvec(i)+gvec(i)*deta;
gvec(i+1)=gvec(i)+hvec(i)*deta;
hvec(i+1)=hvec(i)-1/2*fvec(i)*hvec(i)*deta;
%
end
%
gvec(total/deta)
%
figure(1)
plot(etavec,fvec)
hold on
plot(etavec,gvec,'r--')
plot(etavec,hvec,'g-.')
%