function h = hfun_new( x, y )
%x = p(:,1); y = p(:,2);
% must accept vectorized input


hS = 145; hT = 50; tS = 8; hSensor = 253;

xT1 = hT/2; yT1 = hS;
T1 = [xT1 yT1];     T4 = [-xT1 yT1];
xT2 = xT1; yT2 = yT1 + hT;
T2 = [xT2 yT2];     T3 = [-xT2 yT2];

T = [T1;T2;T3;T4];
nT = size(T,1);

xS1 = 108; yS1 = 0;
S1 = [xS1 yS1];     S8 = [-xS1 yS1];
xS2 = xS1; yS2 = yS1 + hS;
S2 = [xS2 yS2];     S7 = [-xS2 yS2];
xS3 = xT1 + tS; yS3 = yS2;
S3 = [xS3 yS3];     S6 = [-xS3 yS3];
xS4 = xS3; yS4 = yS2 + hT + tS;
S4 = [xS4 yS4];     S5 = [-xS4 yS4];

S = [S1;S2;S3;S4;S5;S6;S7;S8];
nS = size(S,1);


% preallocate and default value
h = inf(size(x));

% transducer (T)
sizeT = 4;
h(abs(x) <= xT1 & y <= yT2 & y >= yT1) = sizeT;

% line charge (C)

d = 1; % tolerance
sizeC = 2;

% % S3 - S2
% h(abs(x) >= S3(1) & y <= S2(2) + d & y >= S2(2)) = sizeC;
%  
% % S3 - S4
% h(abs(x) => S3(1) & abs(x) <= S3(1) + d & y <= S3(2) & y => S4(2)) = sizeC;
%  
% % S5 - S4
% h(abs(x) <= S4(1) & y <= S4(2) + d & y >= S4(2)) = sizeC;


end

