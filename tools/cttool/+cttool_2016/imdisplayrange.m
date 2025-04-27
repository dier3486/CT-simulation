function hp = imdisplayrange(varargin)
%IMDISPLAYRANGE Display Range tool.
%   IMDISPLAYRANGE creates a display range tool in the current figure.
%   The display range tool shows the display range of the intensity image
%   or images in the figure.
%
%   The tool is a uipanel object, positioned in the lower-right corner of
%   the figure, that contains the text string "Display Range:" followed
%   by the display range values for the image. For an indexed, truecolor,
%   or binary image, the display range is not applicable and is set
%   to empty ([]).
%
%   IMDISPLAYRANGE(H) creates a display range tool in the figure specified
%   by the handle H, where H is an image, axes, uipanel, or figure
%   object. Axes, uipanel, or figure objects must contain at
%   least one image object.
%
%   IMDISPLAYRANGE(HPARENT,HIMAGE) creates a tool in HPARENT that shows
%   the display range of HIMAGE. HIMAGE is a handle to an image object or
%   an array of image object handles. HPARENT is handle to the figure or
%   uipanel object that contains the display range tool.
%
%   HPANEL = IMDISPLAYRANGE(...) returns a handle to the display range tool
%   uipanel.
%
%    Examples
%    --------
%
%        imshow('bag.png');
%        imdisplayrange
%
%        dcm = dicomread('CT-MONO2-16-ankle.dcm');
%        subplot(1,2,1), imshow(dcm);
%        subplot(1,2,2), imshow(dcm,[]);
%        imdisplayrange
%
%    See also IMTOOL.

%   Copyright 2004-2010 The MathWorks, Inc.

% validate input and hg objects
[h,parent] = parseInputs(varargin{:});

if strcmp(get(parent,'Type'),'figure')
    parentIsFigure = true;
else
    parentIsFigure = false;
end

imageHandles = imhandles(h);
if isempty(imageHandles)
    error(message('images:imdisplayrange:noImageInFigure'))
end

hFig = ancestor(h,'Figure');
% hFig will be a cell array if imageHandles is an array, even though
% imageHandles may belong to the same figure.
if iscell(hFig)  
    hFig = hFig{1};
end

% initialize for function scope
callbackID = [];

% calculate initial position of tool
units = 'Pixels';
parentPos = getpixelposition(parent);
visibility = 'off';
posPanel = [1 1 parentPos(3) 20];
fudge = 2;

% make sure tool color matches parent
if parentIsFigure
    backgrndColor = get(parent,'Color');
else
    backgrndColor = get(parent,'BackgroundColor');
end

% create panel and label objects
hPanel = uipanel('Parent',parent,...
    'Units',units, ...
    'Visible',visibility,...
    'Tag','displayrange panel',...
    'Bordertype','none',...
    'Position',posPanel,...
    'BackgroundColor', backgrndColor);

labelString = strcat(getString(message('images:imdisplayrangeUIString:labelString')),...
    ':');
hLabel = uicontrol('Parent',hPanel,...
    'Style','text',...
    'String',labelString, ...
    'Tag','displayrange label',...
    'Units','pixels',...
    'Visible',visibility,...
    'BackgroundColor',backgrndColor);

% position label on the left side of panel
labelExtent = get(hLabel,'Extent');
posLabel = [posPanel(1) posPanel(2) labelExtent(3) ...
    labelExtent(4)];
set(hLabel,'Position',posLabel);

% create display range text object
defaultString = ['[',getString(message('images:imdisplayrangeUIString:defaultString')),']'];
hRange = uicontrol('Parent',hPanel,...
    'Style','text',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'BackgroundColor',backgrndColor,...
    'BusyAction','queue',...
    'Enable', 'inactive',...
    'Visible',visibility,...
    'Interruptible','off',...
    'Tag','displayrange text');

% get our image model
imageModels = getimagemodel(imageHandles);

% create actual display range string
if numel(imageModels) == 1
    createStaticDisplayRange;
else
    createDynamicDisplayRange;
    % do this for display purposes
    set(hRange,'String',defaultString);
end

% link visibility of hPanel and its children
hlink = linkprop([hPanel hLabel hRange],...
    'Visible');
setappdata(hPanel,'linkToChildren',hlink);
clear hlink;

% reposition panel and set visible 'on'
setPanelPosition
set(hPanel,'Position', [panelLeft 1 panelWidth panelHeight]);
set(hPanel,'Visible','on');

% Make sure panel will resize with it's figure parent
if isequal(parent,hFig)
    iptaddcallback(hFig,'ResizeFcn',@resizePanel);
    if strcmp(get(parent, 'Visible'),'on')
        figure(parent);
    end
end

% Make sure imdisplayrange can update itself if the image changes
reactToImageChangesInFig(imageHandles,hPanel,@reactDeleteFcn,...
    @reactRefreshFcn);
registerModularToolWithManager(hPanel,imageHandles);

% Set optional output argument
if nargout > 0
    hp = hPanel;
end


    %-------------------------------
    function reactDeleteFcn(obj,evt) %#ok<INUSD>
        
        if ishghandle(hPanel)
            delete(hPanel);
        end
        
    end


    %--------------------------------
    function reactRefreshFcn(obj,evt) %#ok<INUSD>
        
        % remove previous callback
        iptremovecallback(hFig,'WindowButtonMotionFcn',callbackID);
        
        % refresh image models
        imageModels = getimagemodel(imageHandles);

        % create new display range text objects
        if numel(imageModels) == 1
            createStaticDisplayRange;
        else
            createDynamicDisplayRange;
            set(hRange,'String',defaultString);
        end
        
    end


    %--------------------------------
    function createStaticDisplayRange

        % allow function scope
        [formatNumber, containsFloat,needsExponent] = ...
            getNumberFormatFcn(imageModels);

        sethRangeDefaultString;
        sethRangePosition;
        setDisplayRange;

        % update val if axes clim or image cdatamapping changes
        axesHandle = ancestor(imageHandles,'Axes');
        climListener = iptui.iptaddlistener(axesHandle, 'CLim', ...
            'PostSet',@updateVal);

        cdataMappingListener = iptui.iptaddlistener(imageHandles, ...
            'CDataMapping','PostSet',@updateVal);

        % make sure listeners are deleted when tool goes away
        setappdata(hPanel,'updateListeners',...
            [climListener cdataMappingListener]);
        % clear unused references to listeners
        clear climListener cdataMappingListener;


        %------------------------------
        function sethRangeDefaultString

            imageClass = getClassType(imageModels);
            switch imageClass
                case 'double'
                    string = formatNumber(0);
                    if ~containsFloat && ~needsExponent
                        % make string big enough in case imageModel
                        % contains int16-like numbers
                        string = '00000';
                    end
                case {'single','logical'}
                    string = formatNumber(0);
                otherwise
                    % assume integer data type
                    maxValInClass = intmax(imageClass);
                    string = sprintf('%d',maxValInClass);
                    string = repmat('0',size(string));
            end

            defaultString = sprintf('[-%s +%s]',string,string);
            set(hRange,'String',defaultString);
            
        end % sethRangeDefaultString

        %-----------------------
        function setDisplayRange
            
            range = getDisplayRange(imageModels);

            % if the image is an integer type, we round the display range
            % values since this is what MATLAB does internally when
            % rendering images.  Although, the value returned by when
            % getting the CLim might be different, we feel the value
            % displayed accurately represents what is shown and this
            % should be OK.
            if ~containsFloat
                range = round(range);
            end

            if ~isempty(range)
                rangestr = sprintf('[%s %s]',formatNumber(range(1)), ...
                    formatNumber(range(2)));
            else
                rangestr = '[ ]';
            end
            set(hRange,'String',rangestr);
            
        end % setDisplayRange


        %--------------------------------
        function updateVal(obj,eventData) %#ok<INUSD>
            
            setDisplayRange;
            extent = get(hRange,'Extent');
            oldPos = get(hRange,'Position');
            if extent(3) > oldPos(3)
                sethRangePosition;
                setPanelPosition;
            end
            
        end % updateVal
        
    end % createStaticDisplayRange


    %---------------------------------
    function createDynamicDisplayRange

        axesHandles = ancestor(imageHandles,'Axes');
        if iscell(axesHandles)
            axesHandles = [axesHandles{:}]';
        end

        callbackID = iptaddcallback(hFig,'WindowButtonMotionFcn', ...
            @showDisplayRange);
        set(hRange,'DeleteFcn',{@removeCallback,hFig,callbackID});

        % set string of hRange so it will display nicely
        % for common cases.
        set(hRange,'String','[0.00E-000 0.00E-000]');

        sethRangePosition;

        %---------------------------------------
        function showDisplayRange(obj,eventdata) %#ok<INUSD>

            index = findAxesThatTheCursorIsOver(axesHandles);
            if ~isempty(index)
                imModel = imageModels(index);
                rangestr = num2str(getDisplayRange(imModel));
                if ~isempty(rangestr)
                    rangestr = sprintf('[%s]',rangestr);
                else
                    %target image is not an intensity image
                    rangestr = '[ ]';
                end

            else
                rangestr = defaultString;
            end

            set(hRange,'String',rangestr);

        end % showDisplayRange
        
    end % createDynamicDisplayRange

    %-------------------------
    function sethRangePosition

        rangeExtent = get(hRange,'Extent');
        posRange = [posLabel(1)+posLabel(3) posPanel(2) rangeExtent(3) ...
            rangeExtent(4)];
        set(hRange,'Position',posRange);
        
    end

    %------------------------
    function setPanelPosition
        
        parentPos = getpixelposition(parent);
        posRange = get(hRange,'Position');
        panelWidth = posLabel(3) + posRange(3) + fudge;
        panelHeight = max([posLabel(4) posRange(4)]);
        panelLeft = parentPos(3) - panelWidth;
        panelLeft = fixLeftPosIfOnMac(panelLeft);
        setpixelposition(hPanel, [panelLeft 1 panelWidth panelHeight]);
        
    end


    %----------------------------
    function resizePanel(obj,evt) %#ok<INUSD>

        hfigPos = get(hFig,'Position');
        if ishghandle(hPanel)
            panelPos = getpixelposition(hPanel);
            left = hfigPos(3) - panelPos(3);
            left = fixLeftPosIfOnMac(left);
            setpixelposition(hPanel, [left 1 panelPos(3) ...
                panelPos(4)]);
        end
    end

end % imdisplayrange


%-----------------------------------------------------
function removeCallback(obj,eventData,hFig,callbackID) %#ok<INUSL>

if ishghandle(hFig)
    iptremovecallback(hFig,'WindowButtonMotionFcn',callbackID);
end

end


%--------------------------------------
function left = fixLeftPosIfOnMac(left)

% need to move the panel over a little on the mac so that the mac
% resize widget doesn't obstruct view.
if ismac
    left = left - 7;
end

end


%------------------------------------------
function [h,parent] = parseInputs(varargin)

narginchk(0,2);

switch nargin
    case 0
        %IMDISPLAYRANGEPANEL
        h = get(0, 'CurrentFigure');
        if isempty(h)
            error(message('images:imdisplayrange:noCurrentFigure'))
        end
        parent = h;

    case 1
        %IMDISPLAYRANGEPANEL(H)
        h = varargin{1};
        iptcheckhandle(h,{'image','axes','figure'},mfilename,'H',1);
        parent = ancestor(h,'Figure');

    case 2
        %IMDISPLAYRANGEPANEL(HPARENT,HIMAGE)
        parent = varargin{1};
        iptcheckhandle(parent,{'figure','uipanel','uicontainer'},mfilename, ...
            'HPARENT',1);
        h = varargin{2};
        checkImageHandleArray(h,mfilename);
end

end
