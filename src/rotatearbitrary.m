function [xn yn zn] = rotatearbitrary(x,y,z,a,b,c,u,v,w,th)
% Given by Marius Birsanu

L = u^2+v^2+w^2;

xn = (a*(v^2+w^2)+u*(-b*v-c*w+u*x+v*y+w*z)+((x-a)*(v^2+w^2)+u*(b*v+c*w-v*y-w*z))*cos(th)+sqrt(L)*(b*w-c*v-w*y+v*z)*sin(th))/L;
yn = (b*(u^2+w^2)+v*(-a*u-c*w+u*x+v*y+w*z)+((y-b)*(u^2+w^2)+v*(a*u+c*w-u*x-w*z))*cos(th)+sqrt(L)*(-a*w+c*u+w*x-u*z)*sin(th))/L;
zn = (c*(u^2+v^2)+w*(-a*u-b*v+u*x+v*y+w*z)+((z-c)*(u^2+v^2)+w*(a*u+b*v-u*x-v*y))*cos(th)+sqrt(L)*(a*v-b*u-v*x+u*y)*sin(th))/L;