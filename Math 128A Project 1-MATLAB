function [rts] = quarti26644262(C)
%This function figures out the solutions of a quatric function
%INPUT: the coefficients for x^3, x^2, and x written as a vector form
%OUTPUT: the vector of the solutions to f
%
format long
a = C(1); %initial coefficients given
b = C(2);
c = C(3);
f = @(x) x.^4 + a*x.^3 + b*x.^2 + c*x - 1; %Defines the function f 
%
%
%to find the FIRST positive real rt, use bisection
%
a1 = 0;
b1 = 1+abs(a)+abs(b)+abs(c)+1;
TOL = 10^-10;
intervalLength = b1-a1;
while intervalLength > TOL
    c1 = a1 + (b1-a1) / 2;
    fc1 = f(c1);
    if fc1 == 0
        %disp(c1);
        break;
    elseif sign(f(a1))*sign(f(c1)) > 0
        a1 = c1;
    else 
        b1 = c1;
    end
    intervalLength = b1-a1;
end 
r1 = c1; %first positive rt
%
%
%to find the SECOND negative real rt, use bisection
%
a2 = -max(1, 1+abs(a)+abs(b)+abs(c)+1);
b2 = 0;
TOL = 10^-10;
intervalLength = b2-a2;
while intervalLength > TOL
    c2 = a2 + (b2-a2) / 2;
    fc2 = f(c2);
    if fc2 == 0
        %disp(c2);
        break;
    elseif sign(f(a2))*sign(f(c2)) > 0
        a2 = c2;
    else 
        b2 = c2;
    end
    intervalLength = b2-a2;
end 
r2 = c2; %second rt







%to find the THIRD AND FOURTH solution, use Horner's Method 
a = [1  C(1) C(2) C(3) -1];
x0 = r1;
b(1) = a(1);
for i = 2:length(a)
    b(i) = a(i) + x0*b(i-1);
end
%y = b(length(a)) -- ignore the remainder
b = b(1:length(b)-1); %this is the coefficients for the cubic function
%
%
%cubic fuction:
s=b(1);
t=b(2);
u=b(3);
v=b(4);
f2 = @(x) s*x.^3 + t*x.^2 + u* x + v;
%
%Use Horner's Method again to find a quadratic function:
z = [s t u v];
x0 = r2;
b(1) = z(1);
for i = 2:length(z)
    b(i) = z(i) + x0*b(i-1);
end
b = b(1:length(b)-1); %this is the coefficients for the quadratic function
%
%quadratic function: 
s=b(1);
t=b(2);
u=b(3);
f3 = @(x) s*x.^2 + t*x + u; 
d = sqrt(t^2-4*s*u); %use the quadratic formula
if t >0
    r3 = (-t-d)/(2*s);
elseif d == 0
    r3 = r4;
else 
    r3 = (-t+d)/(2*s);
end
r4 = u/(s*r3);
%
%do newton's method on r3 
f1 = @(x) x.^4 + C(1)*x.^3 + C(2)*x.^2 + C(3)*x - 1;
df1 = @(x) ((4*x + 3*C(1))*x + 2*C(2))*x + C(3);
df2 = @(x) (12*x+6*C(1))*x+2*C(2);
p0 = r3;
TOL = 10^-6;
i = 1;
N0 = 10000;
C;
while i <= N0
    p = p0 - f1(p0)*df1(p0)/[(df1(p0))^2-f1(p0)*df2(p0)];
    if (abs(p-p0)<TOL)
        %fprintf('Took %i iterations \n', i);
        break;
    end
    i = i+1;
    p0 = p;
end
newr3 = p;
%
%do newton's method on r4 
f1 = @(x) x.^4 + C(1)*x.^3 + C(2)*x.^2 + C(3)*x - 1;
df1 = @(x) ((4*x + 3*C(1))*x + 2*C(2))*x + C(3);
p0 = r4;
TOL = 10^-6;
i = 1;
N0 = 10000;
C;
while i <= N0
    p = p0 - f1(p0)*df1(p0)/[(df1(p0))^2-f1(p0)*df2(p0)];
    if (abs(p-p0)<TOL)
        %fprintf('Took %i iterations \n', i);
        break;
    end
    i = i+1;
    p0 = p;
end
newr4 = p;
[rts] = [r1; r2; newr3; newr4];



return;
end
