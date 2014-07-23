function gh=ghest_aa(pc,ac,aa,ts,ai)

% GHESTNEW. Calculates GH from regression equations 
% according to Meskers et al (1998),
% but rewritten as a regression equation based on AA, TS and AI as anatomical 
% landmarks

%disp(['Warning in ghest_aa: data must be provided in millimeters!!'])


Rsca=asscap(aa,ts,ai);
Osca=(aa);

pc=Rsca'*(pc-Osca);
ac=Rsca'*(ac-Osca);
aa=Rsca'*(aa-Osca);
ts=Rsca'*(ts-Osca);
ai=Rsca'*(ai-Osca);

laapc=norm(aa-pc);
ltsai=norm(ts-ai);
laipc=norm(ai-pc);
ltspc=norm(ts-pc);

scx=[1 ts(1) laipc]';
scy=[1 ac(2) pc(3)]';
scz=[1 laapc ts(1)]'; 

thx=[26.896   0.614  0.295];
thy=[-16.307   0.825   0.293];
thz=[-1.740   -0.899   -0.229];

GHx = thx*scx;
GHy = thy*scy;
GHz = thz*scz;

gh=[GHx;GHy;GHz]

gh=(Rsca*gh)+Osca;

