% Find the intercept that the true PF intersects each objective axis
function [intercept,pop]=findIntercept(arch,pop,Zmin)
    merged=unique([arch,pop],'rows');
    popObj=merged.objs-repmat(Zmin,length(merged),1);
    [N,M]  = size(popObj);
    % Detect the extreme points and potentially extreme points
    w = 1e-3+eye(M);
    archObj=arch.objs-repmat(Zmin,length(arch),1);
    
    for i = 1 : M
        [~,I]=max(popObj./repmat(w(i,:),N,1),[],2);  
        Potential=find(I==i);
        Y=popObj(Potential,i);
        if (length(Potential)>2)            
            intercept(i)=min(Y);
        
        elseif isempty(Potential)
            
            intercept(i)=max(archObj(:,i),[],1);
        else 
            intercept(i)=max(popObj(:,i),[],1);

        end
    end 
    intercept=intercept+Zmin;
end