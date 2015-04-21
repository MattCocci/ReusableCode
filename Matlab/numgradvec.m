% function [g badg] = numgradvec(fcn,x,varargin)
% modification of numgrad (Sims' program) when the function is a vector and not a scalar - generally used to 
% get second derivatives when the analytical first derivative is available
% note: the output of fcn must be a 1xk (usually 1xn) vector, and the x vector must be nx1

function [g, badg] = numgradvec(fcn,x,varargin)
delta = 1e-6;
n=length(x);
tvec=delta*eye(n);

f0 = eval([fcn '(x,varargin{:})']);
g=zeros(n,size(f0,2));


badg=0;
for i=1:n
   scale=1;  % originally 1
   
   if size(x,1)>size(x,2)
      tvecv=tvec(i,:);
   else
      tvecv=tvec(:,i);
   end
   % g0 is 1xk
   g0 = (eval([fcn '(x+scale*tvecv'', varargin{:})']) - f0) ...
         /(scale*delta);
   if max(abs(g0))< 1e15
      g(i,:)=g0;
   else
      disp('bad gradient ------------------------') % Jinill Kim
      g(i,:)=0;
      badg=1;
   end
end

