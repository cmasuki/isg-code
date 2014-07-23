function [DATA0]=loktrans(DATA,R,v,a)
% program for the transformation of local data to a global coordinate system
% R,v,a are the rotation matrix, translation vector and centerpoint stemming
% from Velpaus.m
% based on the paper by Veldpaus et al (1988)

[m,n]=size(DATA);
vmat = v';
amat = a';
for i = 2:m,
  vmat = [vmat;v'];
  amat = [amat;a'];
end
DATA0 = amat + vmat + [R*(DATA-amat)']';
