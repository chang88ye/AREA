function AREA(Global)
% <algorithm> <A-G>
% AREA: Adaptive reference set guided evolutionary algorithm

%--------------------------------------------------------------------------
% Implemented by Dr Shouyong Jiang, last modified on 15/12/2019
% Please acknowledge the use of this code and reference 
% "Jiang et al. AREA: An adaptive reference-set based evolutionary algorithm 
% for multiobjective optimisation. Information Sciences. Online:
% https://doi.org/10.1016/j.ins.2019.12.011".
%--------------------------------------------------------------------------

%% set default operator
operator=Global.operator;
% operator=@DE;

%% Generate the weight vectors
[W,Global.N] = UniformPoint(Global.N,Global.M);
T = ceil(sqrt(Global.N));
archSize=ceil(1.2*Global.N);

%% Detect the neighbours of each solution
B = pdist2(W,W);
[~,B] = sort(B,2);
B = B(:,1:T);

W0=W; B0=B; W1=W;
%% Generate random population
Population = Global.Initialization();

FrontNo= NDSort(Population.objs,Population.cons,Global.N);
Archive= Population(FrontNo==1);
Z = min(Population.objs,[],1);

%% Optimization
while Global.NotTermination(Archive)
    
    % find intercept
    [Intercept]=findIntercept(Archive, Population,Z);
    
    % update reference set and neighbour
    Freq=floor(0.05*Global.maxgen);
    if mod(Global.gen,Freq)==0
        if mod(Global.gen/Freq,2)==1
            W=W1;
            if length(Archive)>0.5*Global.N
                [Population, W]=updateReferenceSet(Population, Archive, W, Z, Intercept);
                B = pdist2(W,W);
                [~,B] = sort(B,2);
                B = B(:,1:T);
                W1=W;
            end
        else
            W=W0;
            B=B0;
        end
        Population=updatePop(Population, Archive, Intercept, Z, W);
    end
    
    % For each solution  
    [dist,id] = pdist2(Archive.objs,Population.objs,'euclidean','Smallest',1);
    [uid,~,ic]=unique(id);
    dist2=pdist2(Archive.objs,Archive(uid).objs,'euclidean','Smallest',Global.M+1);
    crowding=prod(dist2(2:end,:),1);
    dist=dist+crowding(ic);
    Prob = dist/max(dist)+0.2;
    
    Offspring=INDIVIDUAL();
    for i = randperm(Global.N) %
        if rand<Prob(i)
            P = B(i,randperm(size(B,2)));
        else
            P= randperm(Global.N);
        end
        
        P=P(1:2);
        % Generate an offspring
        if isequal(operator,@DE)
            P=[i, P];
        end
        
        child = Global.Variation(Population(P),1, operator);
        
        % Update the ideal point
        Z = min(Z,child.obj);
        
        % Update solution of best matched reference
        [~,bestR]=min(fitnessMat((child.obj-Z)./(Intercept-Z),W));
        
        P=B(bestR,:);     
        g_old= fitness((Population(P).objs - Z)./(Intercept-Z),W(P,:));
        g_new= fitness(repmat((child.objs-Z)./(Intercept-Z), length(P),1),W(P,:));
        
        Population(P(find(g_old>=g_new,2)))=child;
        Offspring(i)=child;
    end
    
    % keep archive of fixed size
    Offspring(NDSort(Offspring.objs,Offspring.cons,1)~=1)=[];
    Archive=[Archive,Offspring];
    [~,IA]=unique([Archive.objs, Archive.cons],'rows');
    Archive=Archive(IA);
    Archive(NDSort(Archive.objs, Archive.cons, 1)~=1)=[];
    
    if (Global.evaluated>=Global.evaluation)
        archSize=Global.N;
    end
    if size(Archive.objs,1)>archSize
        del=truncate(Archive.objs,size(Archive.objs,1)-archSize);
        Archive(del)=[];
    end
    
end
end



function pop=updatePop(pop, arch, zmax, zmin, W)
[N,M]=size(W);
merged=[pop,arch];
[~,IA]=unique(merged.objs,'rows');
if length(IA)<=N
    return;
end

merged=merged(IA);
mergedObj=merged.objs;

mergedObj=(mergedObj-repmat(zmin, length(merged),1))./repmat(zmax-zmin, length(merged),1);


record=true(N,1);
dist2 = pdist2(W-1.0/M,mergedObj);
k=0;
while (k<N)
    [dist,id]=min(dist2,[],1);
    if isinf(dist)
        pop(record)=merged(randperm(sum(record)));
        break;
    end
    for i=1:N
        if (record(i))
            idi=find(id==i);
            [~,s]=min(dist(idi));
            if ~isempty(s)
                pop(i)=merged(idi(s));
                record(i)=false;
                dist2(:,idi(s))=inf;
                dist2(i,:)=inf;
                k=k+1;
            end
        end
        if (k==N)
            break;
        end
    end
end
end