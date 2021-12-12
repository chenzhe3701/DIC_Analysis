% [misAngleD, misAxis, RodriguersVector] = calculate_misorientation_hcp(euler1_in_degree, euler2_in_degree)
%
% Calculate misorientation. Works for HCP
% Need function derotation
%
% Zhe Chen, 2015-08-04 revised.
% Zhe Chen, 2016-06-20, removed the loop on orientation-2.
% Add the Rodrigues vector representation rho
% Also, comment on 2016-06-20, this is the rotation from orientation-2 to
% to orientation-1, but the reference frame is the intermediate frame.
% i.e., the axis is the axis expressed in the INTERMEDIATE system!!!
%
% updated 2017-03-29.  If the two euler-angles sets are from different settings, put an additional variable, such as 'euler-2 setting-2'
% Otherwise

function [misAngleD, misAxis, rho] = calculate_misorientation_hcp(euler1_in_degree, euler2_in_degree, varargin)
euler1_in_degree = euler1_in_degree + 0.00001;  % in case they are the same

Mextra1 = [1 0 0; 0 1 0; 0 0 1];
Mextra2 = [1 0 0; 0 1 0; 0 0 1];
if ~isempty(varargin)
    if strcmpi(varargin{1},'euler-2 setting-2')
        Mextra2 = euler_to_transformation([0 0 0],[-90, 180, 0], [0 0 0]);
    end
end

misAngleD = 361;
phi1 = euler1_in_degree(1);
Mphi1 = [cosd(phi1) sind(phi1) 0; -sind(phi1) cosd(phi1) 0; 0 0 1];

p1 = euler2_in_degree(1);
Mp1 = [cosd(p1) sind(p1) 0; -sind(p1) cosd(p1) 0; 0 0 1];

ii=1;
for iphi = 0:180:180
    phi = euler1_in_degree(2) + iphi;
    Mphi = [1 0 0 ; 0 cosd(phi) sind(phi); 0 -sind(phi) cosd(phi);];
    
    for iphi2 = 0:60:300
        phi2 = (1-2*iphi/180)*euler1_in_degree(3) + iphi2;
        Mphi2 = [cosd(phi2) sind(phi2) 0; -sind(phi2) cosd(phi2) 0; 0 0 1];
        M1 = Mphi2 * Mphi * Mphi1;
        M1 = M1 * Mextra1;
        R1 = M1';
        
        for ip = 0%:180:180
            p = euler2_in_degree(2) + ip;
            Mp = [1 0 0 ; 0 cosd(p) sind(p); 0 -sind(p) cosd(p);];
            
            for ip2 = 0%:60:300
                p2 = (1-2*ip/180)*euler2_in_degree(3) + ip2;
                Mp2 = [cosd(p2) sind(p2) 0; -sind(p2) cosd(p2) 0; 0 0 1];
                M2 = Mp2 * Mp * Mp1;
                M2 = M2 * Mextra2;
                R2 = M2';
                
                R2to1 = R1*M2;
                [thetad, v] = derotation(R2to1);
                
                %                 % This suggest that the loop for orientation-2 is not necessary
                %                 T(ii) = thetad;
                %                 V(:,ii) = v;
                %                 ii = ii+1;
                
                if thetad < misAngleD
                    misAngleD = thetad;
                    misAxis = v;
                end
            end
        end
    end
end

% T = reshape(T,12,[]);
% V = reshape(V, 12*3, []);
rho = misAxis * tand(misAngleD/2);

end