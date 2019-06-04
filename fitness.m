function fit=fitness(popObj,W)
if size(popObj,1)~=size(W,1)
    error('invalid dimension')
end

M=size(W,2);
offset=1.0/M;

R=W-offset;
fit=max(popObj-R,[],2);
end