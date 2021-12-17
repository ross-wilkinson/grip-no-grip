function animateCrankForces(Pxyz,RFxyz)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

set(gcf,...
    'color', [.7 .7 .8], ...                % Figure, color
    'Menubar', 'none', ...
    'numbertitle', 'off');

set(gca,...
    'Units', 'Normalized', ...     % Axes  
    'OuterPosition', [0, 0, 1, .9], ...     % - position
    'TickLabelinterpreter', 'Latex', ...     % - font
    'XLim', [0 12], 'YLim',[-1 1], 'ZLim',[-1 1]); % - limits

info = annotation('TextBox',... % Info Box
    'Position', [0,.9,1,.1], ...  % - position
    'BackgroundColor', [.2 .2 .2], ...  % - color
    'Fontsize', 18, 'Interpreter', 'Latex');     % - font
                
view(120, 23);  
grid on;

% Crank-Force-Arrows

set(info, 'string', 'Crank-Force Arrows',...
    'color', [.5 .5 .5] + rand(1, 3)/2);
       
xlabel('fore-aft')
ylabel('up-down')
zlabel('medio-lateral')
view(2)
xlim auto
ylim auto

for i = 950:960%1:10:length(Pxyz)/2                          % Loop for each arrow
    r     = .5 + cos(i)*cos(i)/2;               % - 'color sphere'
    g     = .5 + cos(i)*sin(i)/2; 
    b     = .5 + sin(i)/2;
    color = [r, g, b];
    
    farrow(Pxyz(1,i), Pxyz(2,i), Pxyz(3,i), RFxyz(1,i), RFxyz(2,i), RFxyz(3,i), color);   % Arrow Function
    pause(.05);
end
pause(1);

end

