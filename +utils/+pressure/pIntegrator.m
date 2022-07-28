function P = pIntegrator(px,py,Pinit,Mask,Increx,Increy)
% Iterative Method for Pressure Estimation
% Pressure solver for time-resolved 2D2C velocity field
%
% Woii junwei.chen@alumnos.uc3m.es 20210922 v1.0
% Experimental Aerodynamics and Propulsion Lab
% of Univeridad Carlos III de Madrid
%
mask = ones(size(Mask) + 2);
mask(2:end-1,2:end-1) = Mask;
mask1 = circshift(mask, 1, 2).*mask;
mask2 = circshift(mask,-1, 2).*mask;
mask3 = circshift(mask, 1, 1).*mask;
mask4 = circshift(mask,-1, 1).*mask;
%
P = zeros(size(px) + 2);
P(2:end-1, 2:end-1) = Pinit;
P(1,:) = P(2,:);   P(end,:) = P(end-1,:);
P(:,1) = P(:,2);   P(:,end) = P(:,end-1);
[xq, yq] = meshgrid([1.5:size(px,2)-0.5], 1:size(px,1));
px = interp2(px, xq, yq, 'spline');
px = [px(1,:);px;px(end,:)];
px1 = [0*px(:,1:2),px,0*px(:,1)];
px2 = [0*px(:,1),px,0*px(:,1:2)];
[xq, yq] = meshgrid(1:size(py,2), [1.5:size(py,1)-0.5]);
py = interp2(py, xq, yq, 'spline');
py = [py(:,1),py,py(:,end)];
py1 = [0*py(1:2,:);py;0*py(1,:)];
py2 = [0*py(1,:);py;0*py(1:2,:)];
% Lambda = 0.1;
% Lambda = 0.25;
for iCount = 1:100000
    % Adaptive relaxation coefficient
    Lambda = 0.15/sqrt(iCount/50000+1) + 0.1;
    PD =(+Increx*px1 - (P - circshift(P, 1, 2))).*mask1 +...
        (-Increx*px2 - (P - circshift(P,-1, 2))).*mask2 +...
        (+Increy*py1 - (P - circshift(P, 1, 1))).*mask3 +...
        (-Increy*py2 - (P - circshift(P,-1, 1))).*mask4;
    P = P + Lambda.*PD;
    P(1,:) = P(2,:); P(end,:) = P(end-1,:);
    P(:,1) = P(:,2); P(:,end) = P(:,end-1);
    if mean(abs(PD(2:end-1,2:end-1)),'all') < 0.0001

%         iCount
        break;
    end
end
iCount;
P = P(2:end-1,2:end-1);
P = reshape(P,[size(P,1)*size(P,2) 1]);
end