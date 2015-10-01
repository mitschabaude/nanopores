function convergencePlot()
% Shows a nice convergence plot using adaptive FEM output
% n   ... Number of Elements
% est ... Estimator
% err ... Error
load('est_h1.mat')
n = N;
est = err;

rate = 1;
prate = '1';
c = 3.7;
f = @(n) exp(-c*n);

loglog(n,est,'r-o','LineWidth',2, 'MarkerFaceColor','r','MarkerSize',5)
%semilogy(n,est,'r-o','LineWidth',2, 'MarkerFaceColor','r','MarkerSize',5)
hold on
%loglog(n2,err,'b-o','LineWidth',2, 'MarkerFaceColor','b','MarkerSize',5)
%#if (nargin>3)
%end
loglog(n,2*est(1)*n(1)^(rate)*n.^(-rate),'k--','LineWidth',1)
%semilogy(n,1e2*est(2)/f(n(2)) * f(n),'k--','LineWidth',1)
%if (nargin>2)
%    legend('Estimator','Error','O(N^{-3/2})')
%end
legend('Estimator',['O(N^{-',prate,'})'])
xlabel('# Elements')
%legend('Newton error',['O(e^{-c*N})'])
%xlabel('Newton iterations')
hold off

end

% Author: Gregor Mitscha-Eibl
