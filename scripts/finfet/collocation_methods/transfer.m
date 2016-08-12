 function [CC]=transfer(x1,x2,y1,y2,z1,z2,spar)
    
 radius=0.5;
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