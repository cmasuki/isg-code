function [R,v,a,e,q]=veldpaus(lok,glob)

% function [R,v,a,e,q]=veldpaus(lok,glob)
% program for the calculation of rotation matrix R en translation vector v
% following the paper by Veldpaus, Woltring en Dortmans (1988)
% each row is interpreted as one point or landmark
% the equation is:
% glob=(R*(lok'-a*ones(1,n))+(a+v)*ones(1,n))', n = number of points
% R  rotation matrix
% v  translation vector
% a  translation vector
% e  summed squared error
% q  (R*(lok'-a*ones(1,n))+(a+v)*ones(1,n))'-glob
% modified by Leonard Rozendaal to allow estimation on the basis of three landmarks

mess = 'unequal local and global matrices!!!!';
[mlok,nlok] = size(lok);
[mglob,nglob] = size(glob);
if mlok ~= mglob,
   mess
   return
end
a = (sum(lok)/mlok)';
p = (sum(glob)/mlok)';
G = zeros(3,3);
for i = 1:mlok,
   G = G + ((glob(i,:)' - p) * (lok(i,:)' - a)')/mlok;
end
[n,Dmu] = eig(G'*G);
[mu2,k] = sort(diag(Dmu));
mu = sqrt(mu2);
detG = det(G);
t = sign(detG);

% Ga = [adj(G)]';
Ga = [cross(G(:,2),G(:,3)) cross(G(:,3),G(:,1)) cross(G(:,1),G(:,2))]';

% mu(1) < mu(2) < mu(3) !
B1 = mu(3) + mu(2) + t*mu(1);
B2 = mu(3)*mu(2) + t*mu(1)*(mu(3) + mu(2));
C = G'*G + B2*eye(3);
R = (Ga + B1*G)*inv(C);
v = p - a;
amat = a;
for i = 2:mlok,
   amat = [amat a];
end
vmat = v;
for i = 2:mlok,
   vmat = [vmat v];
end
pdak = [(amat + vmat + R*(lok'- amat))]';
e = sqrt(sum((pdak - glob)'.^2)');
%
% condities checken
%
A = zeros(3,3);
for i = 1:mlok,
   A = A + ((lok(i,:)' - a)*(lok(i,:)' - a)')/mlok;
end
[nalpha,Dalpha] = eig(A);
alpha = sort(diag(Dalpha));
cC = (alpha(3) + alpha(2))/(alpha(2) + alpha(1));
if cC > 100,
   cC
   disp('possibly bad estimate')
   %pause
end
fmax = 1/(alpha(2) + alpha(1));
if fmax > 0.1,
   fmax
   disp('possibly bad estimate')
   %pause
end
%q=(R*(lok'-a*ones(1,mlok))+(a+v)*ones(1,mlok))'-glob;
q=(R*(lok'-a*ones(1,mlok))+(a+v)*ones(1,mlok))';

if abs(imag(R))>0
%   disp(['R is imaginary !'])
%   disp(R)
%   disp('rounded to real'),pause(2)
   R=real(R);
end

if abs(imag(v))>0
%   disp(['v is imaginary !'])
%   disp(v)
%   disp('rounded to real')
   v=real(v);
end

if abs(imag(a))>0
%   disp(['a is imaginary !'])
%   disp(a)
%   disp('rounded to real')
   a=real(a);
end

if abs(imag(e))>0
%   disp(['e is imaginary !'])
%   disp(e)
%   disp('rounded to real')
   e=real(e);
end

if abs(imag(q))>0
%   disp(['q is imaginary !'])
%   disp(q)
%   disp('rounded to real')
   q=real(q);
end

