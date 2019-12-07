% plot crystal orientation [] on pole figure
% euler angle

function [] = plot_on_PF(input_dir, eulers_d, sample_direction_of_interest)
%%
% The sample_reference_frame I often used is similar to that of imaging conventin, i.e., x=right, y=down, z=inward.
% For pole figure, we need a plot_reference_frame (sample ref rotate again into 'another sample ref')  
% 
% (1) If sample direction of interest = ND = [-sampleZ], then
% right = plotX = sampleX; up = plotY = -sampleY; out = plotZ = -sampleZ
% euler_sys = from plot to sample by phi_sys(1)@z-->phi_sys(2)@x-->phi_sys(3)@z 
% euler_sys = [0 180 0] .  
% 
% (2) If sample direction of interest = RD = sampleX, then
% right = plotX = -sampleY, up = plotY = -sampleZ, out = plotZ = sampleX 
% euler_sys = [0 90 90]
% 
% (3) If sample direction of interest = TD = sampleY, then
% right = plotX = sampleX, up = plotY = -sampleZ, out = plotZ = sampleY 
% euler_sys = [0 -90 0];

switch sample_direction_of_interest
    case 'ND'
        euler_sys = [0 180 0];
    case 'RD'
        euler_sys = [0 90 90];
    case 'TD'
        euler_sys = [0 -90 0];
    otherwise
        error('has to be ND, RD, or TD');
end

% draw frame
figure;
hold on;
axis square;axis off;
axis([-1.2 1.2 -1.2 1.2]);

t = 0:1:360;
cx = cosd(t);
cy = sind(t);
plot(cx,cy,'-k');

lx = -1:0.01:1;
ly = zeros(size(lx));
plot(lx,ly,'-k');

ly = -1:0.01:1;
lx = zeros(size(lx));
plot(lx,ly,'-k');

% draw pole
colors = parula(size(input_dir,1)+1);
for iL = 1:size(input_dir,1)
    for i = 1:size(eulers_d,1)
        
        Euler(1:3) = eulers_d(i,:);
        
        g = euler_to_transformation(Euler, euler_sys, [0 0 0]);
        pole = input_dir(iL,:) * g;
        
        PN = sign(pole(3));
        if PN==0
            PN = -1;
        end
        
        x = -PN * pole(1);
        y = -PN * pole(2);
        z = -PN * pole(3);  % make z negative.
        X = x/(1-z);    % the coord projected to equitorial plane: X/x = 1/(1-z).
        Y = y/(1-z);
        if z == -1  % for the extreme condition...
            X = -10^(-9);
        end
        
        plot(X,Y,'.','MarkerSize',12,'color',colors(iL,:));
        
%         % plot accord to equal surface area
%         if (x^2+y^2<1-0.8^2)
%             plot(X,Y,'.','MarkerSize',12,'color','r');
%         elseif (x^2+y^2<1-0.6^2)
%             plot(X,Y,'.','MarkerSize',12,'color','y');
%         elseif (x^2+y^2<1-0.4^2)
%             plot(X,Y,'.','MarkerSize',12,'color','g');
%         elseif (x^2+y^2<1-0.2^2)
%             plot(X,Y,'.','MarkerSize',12,'color','b');
%         elseif (x^2+y^2<1-0^2)
%             plot(X,Y,'.','MarkerSize',12,'color','m');
%         end
        
    end
end

switch sample_direction_of_interest
    case 'ND'
        text(1.1, 0,'RD','fontsize',18); 
        text(0, 1.1,'TD','fontsize',18);  
    case 'RD'
        text(1.1, 0,'TD','fontsize',18); 
        text(0, 1.1,'ND','fontsize',18); 
    case 'TD'
        text(1.1, 0,'RD','fontsize',18); 
        text(0, 1.1,'ND','fontsize',18); 
    otherwise
        error('has to be ND, RD, or TD');
end


hold off;


