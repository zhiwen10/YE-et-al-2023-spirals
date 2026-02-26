function [ Aremaped ] = remapColor( A,lold,hold)
 Aremaped = zeros(size(A));
%  lold=min(A);
%  hold=max(A);
 for i=1:length(A)
     newVal = 1 + (A(i)-lold)*(255-1)/(hold-lold);
     Aremaped(i) = round(newVal); 
 end
end