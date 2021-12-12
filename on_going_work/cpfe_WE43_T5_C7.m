addChenFunction
load('T5_#7_EbsdToSemForTraceAnalysis.mat')

phiSys = [90, 180, 0];
% align euler angle to sample reference frame ------------ align.  UMich data is actually setting-1 !!!
[phi1,phi,phi2] = align_euler_to_sample(phi1,phi,phi2,'none', phiSys(1),phiSys(2),phiSys(3)); 

[gPhi1,gPhi,gPhi2] = align_euler_to_sample(gPhi1,gPhi,gPhi2,'none', phiSys(1),phiSys(2),phiSys(3));

eulerAligned = 1;

save('WE43_T5_C7_organized.mat','ID','gID','gPhi','gPhi1','gPhi2','phi','phi1','phi2','eulerAligned');  % record if eulerAligned.