
K=3612;

%% compute for DTI
N = 8;
SSD = 18.191;
AIC_DTI = ComputeAIC(N,K,SSD);
BIC_DTI = ComputeBIC(N,K,SSD);

disp(['DTI; AIC: ' num2str(AIC_DTI) ' BIC: ' num2str(BIC_DTI)]);

%% compute for Ball and Stick
N = 6;
SSD = 15.106;
AIC_BallAndStick = ComputeAIC(N,K,SSD);
BIC_BallAndStick = ComputeBIC(N,K,SSD);

disp(['BallAndStick; AIC: ' num2str(AIC_BallAndStick) ' BIC: ' num2str(BIC_BallAndStick)]);

%% compute for Zeppelin and Stick
N = 7;
SSD = 10.8167;
AIC_ZeppelinAndStick = ComputeAIC(N,K,SSD);
BIC_ZeppelinAndStick = ComputeBIC(N,K,SSD);

disp(['ZeppelinAndStick; AIC: ' num2str(AIC_ZeppelinAndStick) ' BIC: ' num2str(BIC_ZeppelinAndStick)]);

%% compute for Zeppelin and Stick Tortuosity
N = 6;
SSD = 11.6052;
AIC_ZeppelinAndStickTortuosity = ComputeAIC(N,K,SSD);
BIC_ZeppelinAndStickTortuosity = ComputeBIC(N,K,SSD);

disp(['ZeppelinAndStickTortuosity; AIC: ' num2str(AIC_ZeppelinAndStickTortuosity) ' BIC: ' num2str(BIC_ZeppelinAndStickTortuosity)]);



function AIC = ComputeAIC(N,K,SSD)
AIC = 2*N + K*log((1/K)*SSD);
end

function BIC = ComputeBIC(N,K,SSD)
BIC = N*log(K) + K*log((1/K)*SSD);
end