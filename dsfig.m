function h = dsfig(f)
% DSFIG Creates or switches to a figure window.
%  
%  function dsfig(f)
%
% Creates or switches to a figure window. When switching to an 
% existing figure window, prevents that window taking the focus.
% If given a string, finds the figure with that name.
% If none exists, creates a figure with that name.
%    
% given
%  
%   f         name or number of figure window
%
%
    
% Category: PLOTTING









if nargin == 0
    h = figure()
    return
elseif f==0
    h = figure()
    return
end

if isnumeric(f)
    if ishandle(f)
        % if the window exists, just set that to
        % be current 
        set(0,'CurrentFigure',f);
        h = gcf; %kchoutag: is this correct?
        return;
    else
        % sadly this will still grab focus.
        % no way around it if the window doesn't exist.
        h = figure(f);
        return;
    end
else        
    fn = findobj('-regexp', 'name', ['\d*:\s*',f, '$']);
    %  if length(fn)==0
    %    fn = findobj('-regexp', 'name', ['\d*:\s*',f]);
    %  end
    if length(fn)>0
        set(0,'CurrentFigure',fn);
        h = gcf;
        return;
    else
        h = figure();
        set(gcf,'name', sprintf('%d: %s', get(gcf,'Number'), f));
        set(gcf,'numbertitle', 'off');
    end
end
