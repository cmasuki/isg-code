  function[S] =  asscap(AA,TS,AI)
% function[S] =  asscap(AA,TS,AI)                                          %
% calculation of local coordinate system of the scapula S.
% Input :                                                                  %
%        aa,ts & ai in column vectors                                      %
% Output :                                                                 %
%        S = [Xs,Ys,Zs] in column vectors                                   %
%                                                                          %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xs = (AA-TS) / norm(AA-TS);
zs = cross(xs,(AA-AI)); zs = zs/norm(zs);
ys = cross(zs,xs);

S = [xs,ys,zs];

