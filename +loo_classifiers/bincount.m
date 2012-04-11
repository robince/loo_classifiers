function C=bincount(data, nR)

C	= zeros(1,nR);
% diff of find of diff trick for counting number of elements
temp 	= sort(data(:));
dtemp	= diff([temp;max(temp)+1]);
count 	= diff(find([1;dtemp]));
indx 	= temp(dtemp>0);
%Pr(indx+1)= count ./ numel(data);	% probability
C(indx+1) = count;
