% remove K members from population
function Del = truncate(PopObj,K)
% Select part of the solutions by truncation
[N,M]=size(PopObj);

% Calculate distance matrix
Distance = inf(N);
if (M<4) % Multiobjective
    Distance = pdist2(PopObj,PopObj);
else % Many-objective (with SDE)
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
end

% Trucation
Del = false(1,size(PopObj,1));
while sum(Del) < K
    Remain   = find(~Del);
    Temp     = sort(Distance(Remain,Remain),2);
    [~,Rank] = sortrows(Temp);
    Del(Remain(Rank(1))) = true;
end
end
