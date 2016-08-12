clc
clear
% ploter=true;
% 
% LEVEL_MAX = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% 
% NSAMPLES = zeros(size(LEVEL_MAX));
% rule_type = 3; 
%  order=4;
% dim_num =4;  
% level_max = LEVEL_MAX(order);
% level_weight = ones(dim_num,1);
% [sparse_XYPHI_new,sparse_weight_new] = test_sgmga(rule_type,dim_num,level_weight,level_max,false);
% sparse_XYPHI_new = (sparse_XYPHI_new+1)/2;
% 



clc
clear
ploter=true;
radius=0.5;
x0=0;y0=0;z0=0;
% x0=start(1);
% y0=start(2);
% z0=start(3);
lw=20;hb=20;wb=20;
number_atom=8;
 

nx=2; lx=lw/nx;
ny=2; ly=hb/ny;
nz=2; lz=wb/nz;
order=3;
randomm=false;

ploter=true;
LEVEL_MAX = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];

NSAMPLES = zeros(size(LEVEL_MAX));
rule_type = 3; 
 

dim_num =number_atom;  
level_max = LEVEL_MAX(order);
dim=dim_num;
level_weight = ones(dim,1);
[sparse_XYPHI_new,sparse_weight_new] = test_sgmga(rule_type,dim,level_weight,level_max,false);
sparsss = (sparse_XYPHI_new+1)/2 ;
sparss=sparsss;
 
  
nsamples=size(sparss,2); 


center=[];
CENTER=[];

for kk=1:nsamples

    jj=1;
 
    spars=sparss(:,kk);

    for l=1:nx
    
        x1=x0+((l-1)*lx);
    
        x2=x0+(l*lx) ;
    
        for j=1:ny
        
            y1=y0+((j-1)*ly);
        
            y2=y0+(j*ly) ;
        
            for k=1:nz
        
                z1=z0+((k-1)*lz);
        
                z2=z0+(k*lz) ;
      
                spar=spars(jj);
            
                [CC]=transfer(x1,x2,y1,y2,z1,z2,spar);
            
                center=[center;CC];
                
                jj=jj+1;
            
            end    
        end
    end
    
    CENTER=[CENTER,center];
    center=[];
end
data=CENTER(:,8);

% for i=1:5
% 
% 
XX=[];YY=[];ZZ=[];
% data=CENTER(:,i);
  for l=0:number_atom-1
      
      XX=[XX;data(3*l+1)];
      YY=[YY;data(3*l+2)];
      ZZ=[ZZ;data(3*l+3)];
      
  end


createfigure(XX, YY, ZZ)
% 
% 
% end