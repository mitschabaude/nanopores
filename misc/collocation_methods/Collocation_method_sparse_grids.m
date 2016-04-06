function [CENTER,weight] = Collocation_method_sparse_grids( start,order,lw,hb, wb,radius,number_atom )

% INPUTS:
 
% order= order of polynomial expansion
% start: starting coordinate of the FET 
% lw: length of the channel (x-axes)
% hb: height of the channel (y-axes)
% wb: width  of the channel  (z-zxes)
% radius: the radius of atoms. 0.5 is recommended.
% number: number of atoms should be 4, 8,16,32 and 64   
% OUTPUT
% center: the coordinates of center of dopant atoms according to collocation method. 

x0=start(1);
y0=start(2);
z0=start(3);
 
switch number_atom

    case 4 
        nx=1;ny=2;nz=2;
        
    case 8 
  nx=2;ny=2;nz=2;
    
    case 16 
    nx=4;ny=2;nz=2;
    
    case 32
         nx=4;ny=4;nz=2;
         
    case 64 
        nx=4;ny=4;nz=4;
end

lx=lw/nx;ly=hb/ny;lz=wb/nz;
  

randomm=false;

ploter=false;
LEVEL_MAX = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];

NSAMPLES = zeros(size(LEVEL_MAX));
rule_type = 1; 
 

dim =number_atom;  
level_max = LEVEL_MAX(order);
level_weight = ones(dim,1);
[sparse_XYPHI_new,sparse_weight_new] = test_sgmga(rule_type,dim,level_weight,level_max,false);
sparss = (sparse_XYPHI_new+1)/2 ;
weight = sparse_weight_new/2^dim;

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

 
     
    
 
 
    function [CC]=transfer(x1,x2,y1,y2,z1,z2,spar)
    
    x1=x1+radius;
    x2=x2-radius;
    X=x1 + spar*(x2-x1);
    
    y1=y1+radius;
    y2=y2-radius;
    Y=y1 + spar*(y2-y1);
    
    z1=z1+radius;
    z2=z2-radius;
    Z=z1 + spar*(z2-z1);
    
    CC=[X;Y;Z];
       
    end




 


end

