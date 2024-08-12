function h = fig(varargin)
% FIG - Creates a figure with a desired size, no white-space, and several other options.
%
%       All Matlab figure options are accepted. 
%       FIG-specific options of the form FIG('PropertyName',propertyvalue,...) 
%       can be used to modify the default behavior, as follows:
%
%       -'units'    : preferred unit for the width and height of the figure 
%                      e.g. 'inches', 'centimeters', 'pixels', 'points', 'characters', 'normalized' 
%                      Default is 'centimeters'
%
%       -'width'    : width of the figure in units defined by 'units'
%                      Default is 14 centimeters
%                      Note: For IEEE journals, one column wide standard is
%                      8.5cm (3.5in), and two-column width standard is 17cm (7 1/16 in)
%
%       -'height'   : height of the figure in units defined by 'units'
%                      Specifying only one dimension sets the other dimension
%                      to preserve the figure's default aspect ratio. 
%
%       -'font'     : The font name for all the texts on the figure, including labels, title, legend, colorbar, etc.
%                      Default is 'Times New Roman' 
%
%       -'fontsize' : The font size for all the texts on the figure, including labels, title, legend, colorbar, etc.
%                      Default is 14pt
%       
%       -'border'   : Thin white border around the figure (compatible with export_fig -nocrop) 
%                      'on', 'off'
%                      Default is 'off' 
%
%   FIG(H) makes H the current figure. 
%   If figure H does not exist, and H is an integer, a new figure is created with
%   handle H.
%
%   FIG(H,...) applies the properties to the figure H.
%
%   H = FIG(...) returns the handle to the figure created by FIG.
%
width=14;
font='Times New Roman';
fontsize=14; 
units='centimeters';
bgcolor='w';
sborder='off';
flag='';
Pindex=[];
%%%%%%%%%%% process optional arguments
optargin = size(varargin,2);
if optargin>0
% check if a handle is passed in
if isscalar(varargin{1}) && isnumeric(varargin{1})
    flag=[flag '1'];
    i=2;
     if ishghandle(varargin{1})==1
        flag=[flag 'i'];
    end
else
    i=1;
end
% get the property values
while (i <= optargin)
    if (strcmpi(varargin{i}, 'border'))
        if (i >= optargin)
            error('Property value required for: %s', num2str(varargin{i}));
        else
            sborder = varargin{i+1};flag=[flag 'b'];
            i = i + 2;
        end
    elseif (strcmpi(varargin{i}, 'width'))
        if (i >= optargin)
            error('Property value required for: %s', num2str(varargin{i}));
        else
            width = varargin{i+1};flag=[flag 'w'];
            i = i + 2;
        end
    elseif (strcmpi(varargin{i}, 'height'))
        if (i >= optargin)
            error('Property value required for: %s', num2str(varargin{i}));
        else
            height = varargin{i+1};flag=[flag 'h'];
            i = i + 2;
        end
    elseif (strcmpi(varargin{i}, 'font'))
        if (i >= optargin)
            error('Property value required for: %s', num2str(varargin{i}));
        else
            font = varargin{i+1};flag=[flag 'f'];
            i = i + 2;
        end
    elseif (strcmpi(varargin{i}, 'fontsize'))
        if (i >= optargin)
            error('Property value required for: %s', num2str(varargin{i}));
        else
           fontsize = varargin{i+1};flag=[flag 's'];
            i = i + 2;
        end
    elseif (strcmpi(varargin{i}, 'units'))
        if (i >= optargin)
            error('Property value required for: %s', num2str(varargin{i}));
        else
            units = varargin{i+1};flag=[flag 'u'];
            i = i + 2;
        end
    elseif (strcmpi(varargin{i}, 'color'))
        if (i >= optargin)
            error('Property value required for: %s', num2str(varargin{i}));
        else
            bgcolor = varargin{i+1};flag=[flag 'c'];
            i = i + 2;
        end
    else
        %other figure properties
        if (i >= optargin)
            error('A property value is missing.');
        else
        Pindex = [Pindex i i+1];
        i = i + 2;
        end
    end
end
end
% We use try/catch to handle errors
try
% creat a figure with a given (or new) handle
if length(strfind(flag,'1'))==1
    s=varargin{1};
    if ishandle(s)==1
    set(0, 'CurrentFigure', s);
    else 
        figure(s);
    end
else
    s=figure;
end
flag=[flag 's'];
% set other figure properties
if ~isempty(Pindex)
    set(s,varargin{Pindex});
end
% set the background color
set(s, 'color',bgcolor);
% set the font and font size
set(s, 'DefaultTextFontSize', fontsize); 
set(s, 'DefaultAxesFontSize', fontsize); 
set(s, 'DefaultAxesFontName', font);
set(s, 'DefaultTextFontName', font);
%%%%%%%%%%% set the figure size
% set the root unit
old_units=get(0,'Units');
set(0,'Units',units);
% get the screen size
scrsz = get(0,'ScreenSize');
% set the root unit to its default value
set(0,'Units',old_units);
% set the figure unit
set(s,'Units',units);
% get the figure's position
pos = get(s, 'Position');
old_pos=pos;
aspectRatio = pos(3)/pos(4);
% set the width and height of the figure
if length(strfind(flag,'w'))==1 && length(strfind(flag,'h'))==1 
    pos(3)=width;
    pos(4)=height;
elseif isempty(strfind(flag,'h'))
    pos(3)=width;
    pos(4) = width/aspectRatio;
elseif isempty(strfind(flag,'w')) && length(strfind(flag,'h'))==1
    pos(4)=height;
    pos(3)=height*aspectRatio; 
end
% make sure the figure stays in the middle of the screen
diff=old_pos-pos;
 if diff(3)<0
 pos(1)=old_pos(1)+diff(3)/2;
     if pos(1)<0
         pos(1)=0;
     end
 end
 if diff(4)<0
 pos(2)=old_pos(2)+diff(4);
    if pos(2)<0
         pos(2)=0;
     end
 end
 
% warning if the given width (or height) is greater than the screen size
if pos(3)>scrsz(3)
warning(['Maximum width (screen width) is reached! width=' num2str(scrsz(3)) ' ' units]);
end
if pos(4)>scrsz(4)
warning(['Maximum height (screen height) is reached! height=' num2str(scrsz(4)) ' ' units]);
end
% apply the width, height, and position to the figure
set(s, 'Position', pos);
if strcmpi(sborder, 'off')
    set(s,'DefaultAxesLooseInset',[0,0,0,0]);
end
    
% handle errors
catch ME
    if isempty(strfind(flag,'i')) && ~isempty(strfind(flag,'s'))
    close(s);
    end
   error(ME.message)
end
s=figure(s);
% return handle if caller requested it.
  if (nargout > 0)
        h =s;
  end