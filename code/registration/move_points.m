function Po=move_points(P,M,center)

Po=zeros(size(P));
% Move center of rotation to origin (0,0))
P(:,1)=P(:,1)-center(1);
P(:,2)=P(:,2)-center(2);

% Rotate
Po(:,1)=P(:,1)*M(1,1)+P(:,2)*M(1,2);
Po(:,2)=P(:,1)*M(2,1)+P(:,2)*M(2,2);

% Move back
Po(:,1)=Po(:,1)+center(1);
Po(:,2)=Po(:,2)+center(2);

% Translate
Po(:,1)=Po(:,1)+M(3,1);
Po(:,2)=Po(:,2)+M(3,2);

        