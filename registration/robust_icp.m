function [TR, TT, E] = robust_icp(q,p,gfun,wfun,rfun,r0,iterOut,iterIn)

Np = size(p,2);

% Transformed data point cloud
pt = p;

% Allocate vector for objective function in every iteration.
E = zeros(iterOut+1,1); 

% Initialize temporary transform vector and matrix.
T = zeros(2,1);
R = eye(2,2);

% Initialize total transform vector(s) and rotation matric(es).
TT = zeros(2,1, iterOut+1);
TR = repmat(eye(2,2), [1,1, iterOut+1]);

kdOBJ = KDTreeSearcher(transpose(q)); % for kdtree

for k=1:iterOut
       
    % Do matching
    %[match,mindist] = match_bruteForce(q,pt);
    [match mindist] = match_kDtree(q,pt,kdOBJ);
    
    p_idx = true(1, Np);
    q_idx = match;
    
    if k == 1
        if r0 == 0
            r = mad(mindist)/0.6745
        else
            r = r0;
        end
        weights = wfun(mindist,r);
        weights = weights ./ sum(weights);
        E(k) = sum(gfun(mindist,r));
    end
    
    %r = median(mindist) - 10
    r = rfun(r);
    
    % Determine weight vector
    %figure; plot(sort(mindist),'.');
    weights = wfun(mindist,r);
    [R,T,EE] = eq_point(q(:,q_idx),pt(:,p_idx),gfun,wfun,r,iterIn);

    % Add to the total transformation
    TR(:,:,k+1) = R*TR(:,:,k);
    TT(:,:,k+1) = R*TT(:,:,k)+T;

    % Apply last transformation
    pt = TR(:,:,k+1) * p + repmat(TT(:,:,k+1), 1, Np);
    T2D = [TR(1,1,k+1) TR(2,1,k+1) 0; TR(1,2,k+1) TR(2,2,k+1) 0; TT(1,1,k+1) TT(2,1,k+1) 1];
    tform = maketform('affine',T2D);
    
    % Objective function 
    mindist = sqrt(sum((pt-q(:,match)).^2,1));
    weights = wfun(mindist,r);
    weights = weights ./ sum(weights);
    E(k+1) = sum(gfun(mindist,r));
    %disp(['Iter ',num2str(k),': EE = ',num2str(EE(end)),', E = ',num2str(E(k+1))]);
    
end

if true
    TR = TR(:,:,end);
    TT = TT(:,:,end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match mindist] = match_bruteForce(q, p)
    m = size(p,2);
    n = size(q,2);    
    match = zeros(1,m);
    mindist = zeros(1,m);
    for ki=1:m
        d=zeros(1,n);
        for ti=1:2
            d=d+(q(ti,:)-p(ti,ki)).^2;
        end
        [mindist(ki),match(ki)]=min(d);
    end
    
    mindist = sqrt(mindist);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [match mindist] = match_kDtree(~, p, kdOBJ)
	[match mindist] = knnsearch(kdOBJ,transpose(p));
    match = transpose(match);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RR,TT,EE] = eq_point(q,p,gfun,wfun,r,iterIn)

m = size(p,2);
n = size(q,2);

pt = p;
RR = eye(2,2);
TT = zeros(2,1);
EE = zeros(1,iterIn);
E = zeros(1,iterIn);  
for k=1:iterIn
    % Compute weights
    mindist = sqrt(sum((pt-q).^2,1));
    weights = wfun(mindist,r);
    if sum(weights) ~= 0
        % normalize weights
        weights = weights ./ sum(weights);

        % find data centroid and deviations from centroid
        q_bar = q * transpose(weights);
        q_mark = q - repmat(q_bar, 1, n);
        % Apply weights
        q_mark = q_mark .* repmat(weights, 2, 1);

        % find data centroid and deviations from centroid
        p_bar = pt * transpose(weights);
        p_mark = pt - repmat(p_bar, 1, m);
        % Apply weights
        p_mark = p_mark .* repmat(weights, 2, 1);

        N = p_mark*transpose(q_mark); % taking points of q in matched order

        [U,~,V] = svd(N); % singular value decomposition

        R = V*U';%V*diag([1 1 det(U*V')])*transpose(U);
        T = q_bar - R*p_bar;

        RR = R*RR;
        TT = R*TT + T;

        %  Apply transformation
        pt = RR * p + repmat(TT, 1, m);
    end
    mindist_post = sqrt(sum((pt-q).^2,1));
    EE(k) = sum(weights.*(mindist_post.^2)) - sum(weights.*(mindist.^2))/2 + sum(gfun(mindist,r));
    E(k) = sum(gfun(mindist_post,r));
end
