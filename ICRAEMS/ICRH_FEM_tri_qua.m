clc;clear;
%网格参数
xmin=-1.0;xmax=1.0;
ymin=0.0;ymax=1.0;
hpml=0.15;                  %pml层厚度
npml=9;                    %pml层层数
%基本常数
fre=35.e6;
P_nH=0.;
mu0=4*3.1415926535e-7;      %真空磁导率
in_J=10000;
omega0=fre*2*pi;

antDphi=0;                 %天线相位差
fprintf("parameter end\n");
inpath='Input\';
fmesh='icrh2D_5_sp8_modif.msh';
fBx='Bx.txt';
fnp='plasmadensity.txt';

%网格划分
tic
[gcoord,nodes3,nodes4,antnna1,antnna2]=mshload([inpath,fmesh]);
fprintf("grid msh load end\n");
toc
ndof=3;                 %每个节点的自由度数
nnode=length(gcoord);   %总节点数
sdof=nnode*ndof;        %总自由度数
nel3=length(nodes3);    %三角形单元个数
nel4=length(nodes4);    %四边形单元个数
nel=nel3+nel4;          %总单元数

eps(:,:)=fepslion(gcoord,P_nH,omega0,[inpath,fnp],[inpath,fBx]);      %计算介电张量
eps(:,:)=fpml(gcoord,nnode,eps,nodes4,xmax,xmin,ymax,ymin,omega0,hpml,npml);

nt=1;
dphi=2*pi/nt;
for it=1:1
    phi=(it-1)*dphi;
    J3=zeros(nel3,3);
    J3(antnna1(:),1)=0;
    J3(antnna1(:),2)=0;
    J3(antnna1(:),3)=in_J*exp(1i*phi);
    J3(antnna2(:),1)=0;
    J3(antnna2(:),2)=0;
    J3(antnna2(:),3)=in_J*exp(1i*(phi+antDphi));
    J3=-1i*omega0.*mu0.*J3;
    J4=zeros(nel4,3);

    tic
    nnel=3;                 %三角形 单元节点数
    edof=nnel*ndof;         %每个单元的自由度数
    index3=feeldof(nodes3,nel3,nnel,ndof);      %得到单元矩阵元素组装到系统矩阵中的索引
    [k,f]=felp2dt3(gcoord,eps,nodes3,J3);
    [kk1,ff1]=feasmbll(k,f,index3,nel3,edof,sdof);
    fprintf("triangle k->kk end\n");
    toc;

    tic
    nnel=4;                     %四边形 单元节点数
    edof=nnel*ndof;         %每个单元的自由度数
    index4=feeldof(nodes4,nel4,nnel,ndof);
    [k,f]=felp2dt4(gcoord,eps,nodes4,J4);      %设置单元矩阵和单元矢量
    [kk2,ff2]=feasmbll(k,f,index4,nel4,edof,sdof);      %装配系统矩阵和系统矢量
    fprintf("quadr k->kk end\n");
    toc;

    kk=kk1+kk2;
    ff=ff1+ff2;

    %[kk,ff]=feaplyc2(kk,ff,bcdof1,bcval1,sdof);

    tic
    fsol=kk\ff;
    fprintf("solve end\n");
    toc

    fsol=permute(reshape(full(fsol),[3,nnode]),[2,1]);
    Ex(:,1)=fsol(:,1);
    Ey(:,1)=fsol(:,2);
    Ez(:,1)=fsol(:,3);
    nodes=zeros(nel,4);
    nodes(:,:)=NaN;
    nodes(1:nel3,1:3)=nodes3';
    nodes(nel3+1:end,1:4)=nodes4'; 
    cobarmax=1000;
    savefile=sprintf('f%-6.2fMHz-phi(0,%-4.3f)',fre/1e6,antDphi);
    savefile=savefile(savefile~=' ');
    if(~exist(savefile,'file'))
        mkdir(savefile);
    end

    x=linspace(xmin,xmax,201);
    y=linspace(ymin,ymax,101);
    [X,Y]=meshgrid(x,y);
    GridEz=griddata(gcoord(:,1),gcoord(:,2),Ez,X,Y);
    FFTGridEz=fftshift(fft(GridEz(:,:),2001,2),2);
    dcroodk=1/(xmax-xmin)*200/2000;
    croodk=[1000:-1:-1000]*dcroodk*2*pi;
    figure(10)
    plot(croodk,(sum(abs(FFTGridEz(23,:)),1)/1e4).^2/20,'linewidth',2)
    xlim([-30,30])
    xlabel('k_|_| (m^-^1)');
    ylabel('power spectrum (a.u.)');
    set(gca,'FontSize',18,'linewidth',2);
    set(gcf, 'unit', 'centimeters', 'position', [10 5 22 13]);
    print(gcf,[savefile,'\Powerspectrum.jpg'],'-djpeg','-r600');
    hold on


    fig1=figure(1);
    p1=patch('Faces',nodes,'Vertices',gcoord(:,:),'FaceVertexCData',real(Ez),'FaceColor','interp');
    h1=colorbar;caxis([-cobarmax,cobarmax]);set(p1,{'LineStyle'},{'none'}) 
    h1.Label.String='Ez (V/m)';h1.Label.FontSize=20;
    name1=sprintf("%s\\realEz%03d.jpg",savefile,it);set(gca,'fontsize',20);
    xlabel('X (m)');ylabel('Y (m)');
    set(gcf, 'unit', 'centimeters', 'position', [10 5 18 13]);
    set(gca, 'position', [0.15 0.2 0.6 0.7]);
    saveas(fig1,name1);close(fig1);
    fig2=figure(2);
    p2=patch('Faces',nodes,'Vertices',gcoord(:,:),'FaceVertexCData',real(Ey),'FaceColor','interp');
    h1=colorbar;caxis([-cobarmax,cobarmax]);set(p2,{'LineStyle'},{'none'}) 
    name2=sprintf("%s\\realEy%03d.jpg",savefile,it);set(gca,'fontsize',20);
    xlabel('X (m)');ylabel('Y (m)');
    h1.Label.String='Ey (V/m)';h1.Label.FontSize=20;
    set(gcf, 'unit', 'centimeters', 'position', [10 5 18 13]);
    set(gca, 'position', [0.15 0.2 0.6 0.7]);
    saveas(fig2,name2);close(fig2);
    fig3=figure(3);
    p3=patch('Faces',nodes,'Vertices',gcoord(:,:),'FaceVertexCData',real(Ex),'FaceColor','interp');
    h1=colorbar;caxis([-cobarmax,cobarmax]);set(p3,{'LineStyle'},{'none'}) 
    name3=sprintf("%s\\realEx%03d.jpg",savefile,it);set(gca,'fontsize',20);
    xlabel('X (m)');ylabel('Y (m)');
    h1.Label.String='Ex (V/m)';h1.Label.FontSize=20;
    set(gcf, 'unit', 'centimeters', 'position', [10 5 18 13]);
    set(gca, 'position', [0.15 0.2 0.6 0.7]);
    saveas(fig3,name3);close(fig3);
    fig4=figure(4);
    p4=patch('Faces',nodes,'Vertices',gcoord(:,:),'FaceVertexCData',real(eps(:,5)),'FaceColor','interp');
    h1=colorbar;caxis([-1e2,1e2]);set(p4,{'LineStyle'},{'none'}) 
    name4=sprintf("%s\\epszz%03d.jpg",savefile,it);set(gca,'fontsize',20);
    xlabel('X (m)');ylabel('Y (m)');
    h1.Label.String='Eps_z_z';h1.Label.FontSize=20;
    set(gcf, 'unit', 'centimeters', 'position', [10 5 18 13]);
    set(gca, 'position', [0.15 0.2 0.6 0.7]);
    saveas(fig4,name4);close(fig4);
end

function eps=fepslion(gcoord,P_nH,omega0,fnp,fBx)
eps0=8.854187817e-12;       %真空介电常数
%mu0=4*3.1415926535e-7;      %真空磁导率
Vc=2.99792458e8;            %m/s
k0=omega0/Vc;               %波数/m
k02=k0*k0;
charge_e=1.6021892d-19;     %电子电荷量
charge_e2=charge_e*charge_e;
mass_e=9.1095e-31;
mass_iD=2*1836.e0*9.1095e-31;
mass_iH=1836.e0*9.1095e-31;
P_nD=1.-P_nH;
omega02=omega0*omega0;
np1D=load(fnp);
Bx1D=load(fBx);
np=interp1(np1D(:,1),np1D(:,2),gcoord(:,2),'linear','extrap');
Bx=interp1(Bx1D(:,1),Bx1D(:,2),gcoord(:,2),'linear','extrap');
clear np1D Bx1D
omega_e2=charge_e2.*np/(eps0*mass_e);           %电子等离子体频率
omega_iD2=charge_e2.*np*P_nD/(eps0*mass_iD);    %D离子等离子体频率
omega_iH2=charge_e2.*np*P_nH/(eps0*mass_iH);    %H离子等离子体频率
Omega_e2=charge_e2.*Bx.*Bx./(mass_e*mass_e);    %电子回旋频率
Omega_iD2=charge_e2.*Bx.*Bx./(mass_iD*mass_iD);
Omega_iH2=charge_e2.*Bx.*Bx./(mass_iH*mass_iH);
eps_iD=omega_iD2./(omega02-Omega_iD2);
eps_iH=omega_iH2./(omega02-Omega_iH2);
eps_e=omega_e2./(omega02-Omega_e2);
eps(:,1)=(1-omega_e2/omega02-omega_iD2/omega02-omega_iH2/omega02);
eps(:,2)=1-eps_iD-eps_iH-eps_e;
clear omega_e2 omega_iD2 omega_iH2
eps_e2=-(sqrt(Omega_e2)./omega0).*eps_e;
eps_iD2=(sqrt(Omega_iD2)./omega0).*eps_iD;
eps_iH2=(sqrt(Omega_iH2)./omega0).*eps_iH;
eps(:,3)=-1i.*(eps_e2+eps_iD2+eps_iH2);
eps(:,4)=-eps(:,3);
eps(:,5)=eps(:,2);
eps(:,:)=eps*k02;
end

function eps=fpml(gcoord,nnode,eps,nodes4,xmax,xmin,ymax,ymin,omega0,hpml,npml)
% eps0=8.854187817e-12;       %真空介电常数
%mu0=4*3.1415926535e-7;      %真空磁导率
Vc=2.99792458e8;            %m/s
k0=omega0/Vc;               %波数/m
k02=k0*k0;
beta0=1;                    %收缩因子 >=1
pa=1;                       %因子幂次数 1
pb=4;                       %因子幂次数 2-4
ps=4;                       %因子幂次数 2-4
index=sparse([nodes4(1,:),nodes4(2,:),nodes4(3,:),nodes4(4,:)],1,1,nnode,1);
index=find(index);          %pml层内的节点下标
x(:)=gcoord(index(:),1);
y(:)=gcoord(index(:),2);
mashxl(:)= (x-xmin)<=hpml;  %左pml边界层
mashxr(:)= (xmax-x)<=hpml;  %右pml
mashyb(:)= (y-ymin)<=hpml;  %下pml
mashyt(:)= (ymax-y)<=hpml;  %上pml

Vc=2.99792458e8;            %m/s
S=real(eps(index,2))/k02;
D=real(1i*eps(index,3))/k02;
N2R=(S+D)*6;              %平行右旋波折射率平方
N2e=(S.*S-D.*D)./S*6;     %垂直非寻常波折射率平方
alpha0=omega0/2;            %频移因子 >=0
lnR=log(10)*(-(log10(npml-1)-1)/log10(2)-3);       %理论反射率
sigma0=-(ps+1)*Vc*lnR/(2*hpml);   %常数d0 衰减系数 >0

ga=(1-abs(x-xmin)/hpml).^pa.*mashxl+(1-abs(x-xmax)/hpml).^pa.*mashxr+(~(mashxl|mashxr));
gb=(1-abs(x-xmin)/hpml).^pb.*mashxl+(1-abs(x-xmax)/hpml).^pb.*mashxr+(~(mashxl|mashxr));
gs=(1-abs(x-xmin)/hpml).^ps.*mashxl+(1-abs(x-xmax)/hpml).^ps.*mashxr+(~(mashxl|mashxr));
alpha=alpha0*(1-ga);                %频移因子alpha
beta=1+(beta0-1)*gb;               %收缩因子beta
sigma=sigma0*gs;                   %衰减系数sigma
Si=ones(length(index),3);
Si(:,1)=(beta+sigma./sqrt(N2R.')./(alpha+1i*omega0)).*(mashxl|mashxr)+(~(mashxl|mashxr));

ga=(1-abs(y-ymin)/hpml).^pa.*mashyb+(1-abs(y-ymax)/hpml).^pa.*mashyt+(~(mashyb|mashyt));
gb=(1-abs(y-ymin)/hpml).^pb.*mashyb+(1-abs(y-ymax)/hpml).^pb.*mashyt+(~(mashyb|mashyt));
gs=(1-abs(y-ymin)/hpml).^ps.*mashyb+(1-abs(y-ymax)/hpml).^ps.*mashyt+(~(mashyb|mashyt));
alpha=alpha0*(1-ga);                %频移因子alpha
beta=1+(beta0-1)*gb;               %收缩因子beta
sigma=sigma0*gs;                   %衰减系数sigma
Si(:,2)=(beta+sigma./sqrt(N2e.')./(alpha+1i*omega0)).*(mashyb|mashyt)+(~(mashyb|mashyt));

eps(index,1)=eps(index,1).*(Si(:,2).*Si(:,3)./Si(:,1)).*(Si(:,2).*Si(:,3)./Si(:,1));
eps(index,2)=eps(index,2).*(Si(:,1).*Si(:,3)./Si(:,2)).*(Si(:,1).*Si(:,3)./Si(:,2));
eps(index,3)=eps(index,3).*Si(:,1).*(Si(:,1).*Si(:,3)./Si(:,2));
eps(index,4)=eps(index,4).*Si(:,1).*(Si(:,1).*Si(:,2)./Si(:,3));
eps(index,5)=eps(index,5).*(Si(:,1).*Si(:,2)./Si(:,3)).*(Si(:,1).*Si(:,2)./Si(:,3));
end

% function [kk,ff]=feaplyc2(kk,ff,bcdof1,bcval1,sdof)
% ibcdof=sparse(1:sdof,1:sdof,1,sdof,sdof);
% iibcdof1=sparse(bcdof1,bcdof1,1,sdof,sdof);
% iiibcdof=ibcdof-iibcdof1;
% kk=iiibcdof*kk;
% kk=kk+iibcdof1;
% ff(bcdof1(:))=bcval1(:);
% end

function [kk,ff]=feasmbll(k,f,index,nel,edof,sdof)
meshx=reshape(kron(ones(edof,1),index),[edof,edof,nel]);
meshy=permute(meshx,[2,1,3]);
% mesh=reshape(kron(ones(1,edof),index),[edof,nel,edof]);
% meshx=permute(mesh,[1,3,2]);
% meshy=permute(mesh,[3,1,2]);
row=reshape(meshx,[edof*edof*nel,1]);
col=reshape(meshy,[edof*edof*nel,1]);
val=reshape(k,[edof*edof*nel,1]);
kk=sparse(row,col,val,sdof,sdof);
frow=reshape(index,[edof*nel,1]);
fval=reshape(f,[edof*nel,1]);
ff=sparse(frow,1,fval,sdof,1);
end

function [k,f]=felp2dt3(gcoord,eps,nodes3,J)
x1(:)=gcoord(nodes3(1,:),1);y1(:)=gcoord(nodes3(1,:),2);
x2(:)=gcoord(nodes3(2,:),1);y2(:)=gcoord(nodes3(2,:),2);
x3(:)=gcoord(nodes3(3,:),1);y3(:)=gcoord(nodes3(3,:),2);
P1(:)=eps(nodes3(1,:),1);S1(:)=eps(nodes3(1,:),2);D1(:)=eps(nodes3(1,:),3);
P2(:)=eps(nodes3(2,:),1);S2(:)=eps(nodes3(2,:),2);D2(:)=eps(nodes3(2,:),3);
P3(:)=eps(nodes3(3,:),1);S3(:)=eps(nodes3(3,:),2);D3(:)=eps(nodes3(3,:),3);

A(:)=0.5*(x2.*y3+x1.*y2+x3.*y1-x2.*y1-x3.*y2-x1.*y3);
nel=length(A);
AA(1,1,:)=A(:);
Jx(1,1,:)=J(:,1);
Jy(1,1,:)=J(:,2);
Jz(1,1,:)=J(:,3);
x21=x2-x1;x32=x3-x2;x13=x1-x3;
y12=y1-y2;y23=y2-y3;y31=y3-y1;
x21=x21';x32=x32';x13=x13';y12=y12';y23=y23';y31=y31';
PyH2=[x32.*x32,     x32.*x13,   x32.*x21;
      x32.*x13,     x13.*x13,   x13.*x21;
      x32.*x21,     x13.*x21,   x21.*x21];
PyH2=permute(reshape(PyH2,[nel,3,3]),[2,3,1])./(4*AA);
PxH2=[y23.*y23,     y23.*y31,   y23.*y12;
      y23.*y31,     y31.*y31,   y31.*y12;
      y23.*y12,     y31.*y12,   y12.*y12];
PxH2=permute(reshape(PxH2,[nel,3,3]),[2,3,1])./(4*AA);
PxPy=[y23.*x32,     y23.*x13,   y23.*x21;
      y31.*x32,     y31.*x13,   y31.*x21;
      y12.*x32,     y12.*x13,   y12.*x21];
PxPy=permute(reshape(PxPy,[nel,3,3]),[2,3,1])./(4*AA);
HxxH=[2*P1.',   P2.',       P3.';
      P1.',     2*P2.',     P3.';
      P1.',     P2.',       2*P3.'];
HxxH=permute(reshape(HxxH,[nel,3,3]),[2,3,1]).*AA/12;
HyyH=[2*S1.',   S2.',       S3.';
      S1.',     2*S2.',     S3.';
      S1.',     S2.',       2*S3.'];
HyyH=permute(reshape(HyyH,[nel,3,3]),[2,3,1]).*AA/12;
HzzH=HyyH;
HyzH=[2*D1.',   D2.',       D3.';
      D1.',     2*D2.',     D3.';
      D1.',     D2.',       2*D3.'];
HyzH=permute(reshape(HyzH,[nel,3,3]),[2,3,1]).*AA/12;
HzyH=-HyzH;

k(1:3,1:3,:)=PyH2-HxxH;
k(1:3,4:6,:)=-permute(PxPy,[2,1,3]);
k(1:3,7:9,:)=0;
k(4:6,1:3,:)=-PxPy;
k(4:6,4:6,:)=PxH2-HyyH;
k(4:6,7:9,:)=-HyzH;
k(7:9,1:3,:)=0;
k(7:9,4:6,:)=-HzyH;
k(7:9,7:9,:)=PxH2+PyH2-HzzH;
f(1:3,1,:)=Jx.*ones(3,1).*AA./3;
f(4:6,1,:)=Jy.*ones(3,1).*AA./3;
f(7:9,1,:)=Jz.*ones(3,1).*AA./3;
end

function [k,f]=felp2dt4(gcoord,eps,nodes4,J)
x1(:)=gcoord(nodes4(1,:),1);y1(:)=gcoord(nodes4(1,:),2);      %单元内节点坐标
x2(:)=gcoord(nodes4(2,:),1);y2(:)=gcoord(nodes4(2,:),2);
x3(:)=gcoord(nodes4(3,:),1);y3(:)=gcoord(nodes4(3,:),2);
x4(:)=gcoord(nodes4(4,:),1);y4(:)=gcoord(nodes4(4,:),2);
%两点高斯积分
gauss=[-0.57735026918962,0.57735026918962];         %高斯积分样本点
weight=[1.0,1.0];                                   %高斯积分权重
nel=length(x1);
PyN2=zeros(4,4,nel);                %型函数y方向偏导^2积分
PxN2=zeros(4,4,nel);                %x方向偏导^2积分
PxPy=zeros(4,4,nel);                %x方向偏导*y方向偏导积分
NxxN=zeros(4,4,nel);
NyyN=zeros(4,4,nel);
NyzN=zeros(4,4,nel);
NzyN=zeros(4,4,nel);
NzzN=zeros(4,4,nel);
NPxxN=zeros(4,4,nel);
NPxyN=zeros(4,4,nel);
Jx(1,1,:)=J(:,1);
Jy(1,1,:)=J(:,2);
Jz(1,1,:)=J(:,3);
fNx=zeros(1,4,nel);
fNy=zeros(1,4,nel);
fNz=zeros(1,4,nel);
epsxx=reshape(eps(nodes4(:,:),1),[1,4,nel]);
epsyy=reshape(eps(nodes4(:,:),2),[1,4,nel]);
epsyz=reshape(eps(nodes4(:,:),3),[1,4,nel]);
epszy=reshape(eps(nodes4(:,:),4),[1,4,nel]);
epszz=reshape(eps(nodes4(:,:),5),[1,4,nel]);
%二维高斯点
for j=1:length(gauss)
    for i=1:length(gauss)
        s=gauss(i);t=gauss(j);      %高斯积分样本点坐标
        w=weight(i).*weight(j);     %高斯积分权重
        J11(1,1,:)=(-(1-t)*x1+(1-t)*x2+(1+t)*x3-(1+t)*x4)/4;    %雅克比矩阵元素
        J12(1,1,:)=(-(1-t)*y1+(1-t)*y2+(1+t)*y3-(1+t)*y4)/4;
        J21(1,1,:)=(-(1-s)*x1-(1+s)*x2+(1+s)*x3+(1-s)*x4)/4;
        J22(1,1,:)=(-(1-s)*y1-(1+s)*y2+(1+s)*y3+(1-s)*y4)/4;
        detJ(1,1,:)=J11.*J22-J21.*J12;                                  %雅克比矩阵行列式
        N=[(1-s)*(1-t)/4,(1+s)*(1-t)/4,(1+s)*(1+t)/4,(1-s)*(1+t)/4];    %型函数
        PsN=[(t-1)/4,(1-t)/4,(1+t)/4,-(1+t)/4];                         %型函数对s偏导
        PtN=[(s-1)/4,-(1+s)/4,(1+s)/4,(1-s)/4];                         %型函数对t偏导
        PstN=[1,-1,1,-1]/4;

        PyN=-J21.*PsN+J11.*PtN;
        PyN2=PyN2+permute(PyN,[2,1,3]).*PyN./detJ.*w;        %二维高斯积分累加得到 型函数y方向偏导^2积分
        
        PxN=J22.*PsN-J12.*PtN;
        PxN2=PxN2+permute(PxN,[2,1,3]).*PxN./detJ.*w;
        
        PxPy=PxPy+permute(PxN,[2,1,3]).*PyN./detJ.*w;
        
        NT=N.';
        
        NxxN=NxxN+NT.*epsxx.*N.*detJ.*w;
        NyyN=NyyN+NT.*epsyy.*N.*detJ.*w;
        NyzN=NyzN+NT.*epsyz.*N.*detJ.*w;
        NzyN=NzyN+NT.*epszy.*N.*detJ.*w;
        NzzN=NzzN+NT.*epszz.*N.*detJ.*w;
        
        fNx=fNx+Jx.*N.*detJ.*w;
        fNy=fNy+Jy.*N.*detJ.*w;
        fNz=fNz+Jz.*N.*detJ.*w;
        
        PxxN=2*PstN.*J22.*(-J12)./detJ./detJ;
        NPxxN=NPxxN+NT.*PxxN.*detJ.*w;
        PxyN=PstN.*(J21.*J12+J22.*J11)./detJ./detJ;
        NPxyN=NPxyN+NT.*PxyN.*detJ.*w;
    end
end
ans1=NPxxN+PxN2;
ans2=NPxyN+PxPy;
k(1:4,1:4,:)=PyN2-NxxN;
k(1:4,5:8,:)=-PxPy;
k(1:4,9:12,:)=0;
k(5:8,1:4,:)=-permute(PxPy,[2,1,3]);
k(5:8,5:8,:)=PxN2-NyyN;
k(5:8,9:12,:)=-NyzN;
k(9:12,1:4,:)=0;
k(9:12,5:8,:)=-NzyN;
k(9:12,9:12,:)=PxN2+PyN2-NzzN;
f(1:4,1,:)=fNx;
f(5:8,1,:)=fNy;
f(9:12,1,:)=fNz;
end

function index=feeldof(nodes,nel,nnel,ndof)
index=zeros(nnel*ndof,nel);
k=0;
for j=1:ndof
    for i=1:nnel
        start(:)=(nodes(i,:)-1)*ndof;
        k=k+1;
        index(k,:)=start(:)+j;
    end
end
end

