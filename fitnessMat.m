function fit=fitnessMat(popObj,W)
M=size(W,2);
offset=1.0/M;

R=W-offset;

for i=1:size(popObj,1)
    fit(i,:)=max(popObj(i,:)-R,[],2)';
end
end