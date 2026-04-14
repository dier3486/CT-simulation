function hout = imcttool(varargin)
% CTTOOL Display image by imtool.  
% modified from IMTOOL by Dier Zhang
%   IMTOOL opens a new Image Tool in an empty state. Use the File menu options
%   "Open..." or "Import From Workspace..." to choose an image for display.
%
%   IMTOOL(I) displays the grayscale image I.
%
%   IMTOOL(I,[LOW HIGH]) displays the grayscale image I, specifying the display
%   range for I in [LOW HIGH]. The value LOW (and any value less than LOW)
%   displays as black, the value HIGH (and any value greater than HIGH) displays
%   as white. Values in between are displayed as intermediate shades of gray,
%   using the default number of gray levels. If you use an empty matrix ([]) for
%   [LOW HIGH], IMTOOL uses [min(I(:)) max(I(:))]; the minimum value in I
%   displays as black, and the maximum value displays as white.
%
%   IMTOOL(RGB) displays the truecolor image RGB.
%
%   IMTOOL(BW) displays the binary image BW. Values of 0 display as black, and
%   values of 1 display as white.
%
%   IMTOOL(X,MAP) displays the indexed image X with colormap MAP.
%
%   IMTOOL(FILENAME) displays the image contained in the graphics file FILENAME.
%   The file must contain an image that can be read by IMREAD or DICOMREAD or a
%   reduced resolution dataset (R-Set) created by RSETWRITE. If the file
%   contains multiple images, the first one will be displayed. The file must
%   be in the current directory or on the MATLAB path.
%
%   HFIGURE = IMTOOL(...) returns a handle HFIGURE to the figure created by
%   IMTOOL. CLOSE(HFIGURE) closes the Image Tool.
%
%   IMTOOL CLOSE ALL closes all instances of the Image Tool.
%
%   IMTOOL(...,PARAM1,VAL1,PARAM2,VAL2,...) displays the image, specifying
%   parameters and corresponding values that control various aspects of the
%   image display. Parameter names can be abbreviated, and case does not matter.
%
%   Parameters include:
%
%   'Colormap'               2-D, real, M-by-3 matrix specifying the figure
%                            colormap property. Use this parameter to view
%                            grayscale images in false color.
%
%   'DisplayRange'           Two-element vector [LOW HIGH] that controls 
%                            the display range of a grayscale image. See 
%                            above for more details about how to set [LOW HIGH].
%
%                            Including the parameter name is optional, except
%                            when the image is specified by a filename. 
%                            The syntax IMTOOL(I,[LOW HIGH]) is equivalent to
%                            IMTOOL(I,'DisplayRange',[LOW HIGH]).
%                            The parameter name must be specified when 
%                            using IMTOOL with a filename, as in the syntax
%                            IMTOOL(FILENAME,'DisplayRange'[LOW HIGH]).
%
%   'InitialMagnification'   One of two text strings: 'adaptive' or 'fit',
%                            or a numeric scalar value that specifies the 
%                            initial magnification used to display the image. 
%
%                            When set to 'adaptive' then on initial display, 
%                            the entire image is visible. If the image is too
%                            large to display on the screen, IMTOOL displays
%                            the image at the largest magnification that fits 
%                            on the screen.
%
%                            When set to 'fit' IMTOOL scales the entire image
%                            to fit in the window.
%
%                            When set to 100, the image is displayed at 100%
%                            magnification. Similarly, for any other numeric
%                            scalar, so when set to 50, the image is displayed
%                            at 50% magnification. 
%
%                            By default, the initial magnification is set to
%                            the value returned by
%                            IPTGETPREF('ImtoolInitialMagnification').
%
%                            When the image aspect ratio is such that less
%                            than one pixel would be displayed at the
%                            requested magnification in either dimension, a
%                            warning is issued and the image is displayed
%                            at 100%.
%
%   Class Support
%   -------------  
%   A truecolor image can be uint8, uint16, single, or double. An indexed image
%   can be logical, uint8, single, or double. A grayscale image can be uint8,
%   uint16, int16, single, or double. A binary image is of class logical.
%
%   If your grayscale image is single or double, the default display range is
%   [0 1]. If your image's data range is much larger or smaller than the default
%   display range, you may need to experiment with setting the display range to
%   see features in the image that would not be visible using the default
%   display range. For all grayscale images having integer types, the default
%   display range is [intmin(class(I)) intmax(class(I))].
% 
%   Large Data Support
%   ------------------
%   To view very large TIFF images that will not fit into memory you can use
%   RSETWRITE to create a reduced resolution dataset (R-Set) viewable in
%   IMTOOL. R-Sets also improve performance of IMTOOL for large images.
%
%   Related Toolbox Preferences
%   ---------------------------
%   You can use the Image Processing Preferences dialog to set toolbox
%   preferences that modify the behavior of IMTOOL. To access the dialog,
%   select File>Preferences in the MATLAB desktop or Image Tool menu.
%   Additionally you can set these preferences programmatically using the
%   IPTSETPREF function. IMTOOL preferences include:
%
%   - 'ImtoolInitialMagnification' controls the initial magnification for
%     image display. To override this toolbox preference, specify the
%     'InitialMagnification' parameter when you call IMTOOL, as follows:
%        imtool(...,'InitialMagnification',initial_mag).
%
%   - 'ImtoolStartWithOverview' controls whether the Overview tool opens
%     automatically when you open an image using the Image Tool. Possible
%     values:
%        true    - Overview tool opens when you open an image.
%        {false} - Overview tool does not open when you open an image. This
%                  is the default behavior.
%
%   For more information about these preferences, see the reference entry for
%   IPTPREFS.
%
%   Remarks
%   -------
%   IMSHOW is the toolbox's fundamental image display function, optimizing 
%   figure, axes, and image object property settings for image display. IMTOOL
%   provides all the image display capabilities of IMSHOW but also provides 
%   access to several other tools for navigating and exploring images, such as
%   the Pixel Region tool, Image Information tool, and the Adjust Contrast 
%   tool. IMTOOL presents an integrated environment for displaying images and
%   performing some common image processing tasks.
%
%   Examples
%   --------
%       % Display an image from a file
%       imtool('board.tif')  
%
%       % Display an indexed image
%       [X,map] = imread('trees.tif');
%       imtool(X,map) 
%
%       % Display a grayscale image 
%       I = imread('cameraman.tif');
%       imtool(I) 
%
%       % Display a grayscale image, adjust the display range
%       h = imtool(I,[0 80]);
%       close(h)
%
%   See also IMAGEINFO, IMCONTRAST, IMOVERVIEW, IMPIXELREGION, IMREAD,
%            IMSHOW, IPTPREFS, RSETWRITE.

%   Copyright 2004-2015 The MathWorks, Inc.  
%   Copy free modified by Dier Zhang

% Non-default XData and YData are not currently supported

%   'XData'                  Two-element vector X establishes a
%                            nondefault spatial coordinate system, by
%                            specifying the image XData. Note that X can
%                            have more than 2 elements, but only the first
%                            and last elements are actually used.
%
%   'YData'                  Two-element vector Y establishes a
%                            nondefault spatial coordinate system, by
%                            specifying the image YData. Note that Y can
%                            have more than 2 elements, but only the first
%                            and last elements are actually used.

  if ~images.internal.isFigureAvailable()
      error(message('images:imtool:needJavaFigure'))
  end

  % short circuit for IMTOOL CLOSE ALL
  close_all =  nargin == 2 && ischar(varargin{1}) && ischar(varargin{2}) && ...
      strcmpi(varargin{1},'close') && ...
      strcmpi(varargin{2},'all');

  % short circuit for IMTOOL CLOSE ALL
  if close_all
     closeAll
     return
  end
  
  %Create invisible figure
  hFig = figure('Toolbar','none',...
                'Menubar','none',...
                'HandleVisibility','callback',...
                'IntegerHandle','off',...
                'NumberTitle','off',...
                'Tag','imtool',...
                'Visible','off',...
                'WindowStyle',...
                get(0,'FactoryFigureWindowStyle'),...
                'DeleteFcn',@deleteTool);
            
 suppressPlotTools(hFig);

  % Set default 'HitTest','off' for figure children.  This guarantees that
  % our image ButtonDownFcn will not be intercepted unexpectedly (g432132)
  turnOffDefaultHitTestFigChildren(hFig);
 
  % initialize for function scope
  hIm = [];
  imgModel = [];
  showDisplayRange = [];
  image_name = '';
  sp_api = [];
  cb_api = [];

  %Tools
  hOverviewFig = []; 
  hImageInfoFig = [];
  hColormapSelectFig = [];
  hBottomPanel = [];
  hPixelRegionFig = [];
  hScrollPanel = [];
  hRangePanel = [];
  hPixInfoPanel = [];
  hContrastFig = [];  
  
  imtoolModeManager = makeUIModeManager(@makeDefaultModeCurrent);
  [zoomInItem, zoomOutItem, panItem, winLevelItem,cropItem,distanceItem] = deal([]);
  [zoomInTool, zoomOutTool, panTool, winLevelTool, cropTool, distanceTool] = deal([]);
  
  %Keeps list of tools to enable
  contrastItems = [];
  colormapItem = [];
  emptyToolInactiveList = [];
 
  % Needed by printToFig, showImageInfo, etc. 
  [cdata, cdatamapping, clim, map, xdata, ydata,...
   initial_mag,filename] = deal([]);

  % Create menus and toolbar 
  toolbar = createMenusAndToolbar();
  
  minFigWidth = getToolbarWidth(toolbar);
  minFigHeight = 128; % don't scrunch too much
  lastFigPos = []; % used in the ResizeFcn
  imageWidth = [];
  imageHeight = [];
  hImagePage = [];
  
  %Stores a handle to an RSet object if one is passed to imtool
  rset = [];
  
  % Set up modes such that tool and menu items stay in sync
  % Must be defined after creating toolbar and menus
  imtoolModeManager.addMode(zoomInTool, zoomInItem, @makeZoomInModeCurrent, @reactToModeChange);
  imtoolModeManager.addMode(zoomOutTool,zoomOutItem,@makeZoomOutModeCurrent, @reactToModeChange);
  imtoolModeManager.addMode(panTool,    panItem,    @makePanModeCurrent, @reactToModeChange);
  imtoolModeManager.addMode(distanceTool, distanceItem, @makeDistanceModeCurrent, @reactToModeChange);
  imtoolModeManager.addMode(cropTool,   cropItem,   @makeCropModeCurrent);

  id_stream = idStreamFactory('ImtoolInstance');
  tool_number = id_stream.nextId(); 

  % Figure out variable name of image for use in naming of window in
  % nested function addImageToImtool().
  input_image_name = '';
  if nargin >= 1
      input_image_name = inputname(1);
  end
  
  is_imtool_empty = (nargin == 0 || isempty(varargin{1}));
       
  if ~is_imtool_empty
      try
          addImageToImtool(varargin{:});
      catch ME
          % If we failed to open a large TIF image and the user elected not
          % to create an rset, assign hout and return without error.
          if strcmp(ME.identifier,'images:getImageFromFile:OutOfMemTif')
              if nargout > 0
                  close(hFig)
                  hout = [];
              end
              return;
          else
              rethrow(ME)
          end
      end
      
      % Turn on tools by default according to preferences
      if iptgetpref('ImtoolStartWithOverview')
          showOverview
      end
      
      imtoolModeManager.activateDefaultMode()
      
  else
      set(emptyToolInactiveList,'Enable','off');
%       set(hFig,'Name',getString(message('images:imtoolUIString:toolNameNew',tool_number)) );
      set(hFig,'Name',sprintf('CT Tool %d - (new)', tool_number));
      set(hFig,'Visible','on');
  end

  % mouse wheel for images page up/down
  set(hFig, 'windowscrollWheelFcn', @imageScroll);
  
  % Install pointer manager.
  iptPointerManager(hFig);

  figure(hFig) % Bring main tool to front

  if (nargout > 0)
    % Only return handle if caller requested it.
    hout = hFig;
  end
    
   %----------------------------------
   function addImageToImtool(varargin)
    
     try
       % Opens the image
       specificArgNames = {}; % No specific args needed
       
       % parse all arguments
       common_args = images.internal.imageDisplayParsePVPairs(specificArgNames,varargin{:});
       
       filename_specified = ~isempty(common_args.Filename);
       if filename_specified 
           common_args = parseFilename(common_args);
       end
       if ~iscell(common_args.CData)
           cdata = common_args.CData;
       else
           cdata = cell2mat(reshape(common_args.CData, 1, 1, length(common_args.CData)));
       end
       hImagePage = iptui.imagePage(size(cdata, 3));
%        pageindex = hImagePage.pageindex;
       common_args.CData = real(squeeze(cdata(:, :, hImagePage.pageindex, :)));
       common_args = images.internal.imageDisplayValidateParams(common_args);
       
       cdatamapping = common_args.CDataMapping;
       clim = common_args.DisplayRange;
       initial_mag = common_args.InitialMagnification;
       map = common_args.Map;
       
       filename = common_args.Filename;
       xdata = common_args.XData;
       ydata = common_args.YData;
       
       % imageDisplayParseInputs is more permissive regarding data types
       % than what imtool and its uitools allow.  Throw an error if cdata is
       % int8, int32, or uint32.
       classCdata = class(cdata);
       if any(strncmp(classCdata, {'int8','uint32','int32'}, ...
               length(classCdata)))
           error(message('images:imtool:invalidType', classCdata))
       end
       
       if isempty(initial_mag)
         initial_mag = iptgetpref('ImtoolInitialMagnification');
       else
         initial_mag = images.internal.checkInitialMagnification(initial_mag,{'fit','adaptive'},...
                                                 mfilename,'INITIAL_MAG', ...
                                                 []);
       end
       
       if isempty(image_name)
         image_name = getImageName(filename, input_image_name);
       end
%        imtool_name = getString(message('images:imtoolUIString:toolNameWithImageName',tool_number,image_name));
       imtool_name = sprintf('CT Tool %d - %s', tool_number,image_name);
       
       set(hFig,'Name', imtool_name );

       % Check if XData or YData are non-default, warn and reset
       imageWidth  = size(cdata,2);
       imageHeight = size(cdata,1); 
       defaultXData = [1 imageWidth];
       defaultYData = [1 imageHeight];   
       
       isXDataDefault = isequal(xdata,defaultXData);
       isYDataDefault = isequal(ydata,defaultYData);
       
       % do not enforce xdata/ydata constraint for rsets
       if isempty(rset) && (~isXDataDefault || ~isYDataDefault) && ...
               (~isempty(xdata) || ~isempty(ydata))
           
           if ~isXDataDefault
               xdata = [1 imageWidth];
               msgXData = [getString(message('images:imtool:specifiedNonDefaultData',...
                   'XData',1,imageWidth)) '\n'];
           else
               msgXData = '';
           end
           
           if ~isYDataDefault
               ydata = [1 imageHeight];
               msgYData = getString(message('images:imtool:specifiedNonDefaultData',...
                   'YData',1,imageHeight));
           else
               msgYData = '';
           end
           
           warning(message('images:imtool:nonDefaultXDataOrYData', msgXData, msgYData))
       end
       
       hAx = axes('Parent', hFig);
       
       warnstate = warning('off','images:initSize:adjustingMag');
       
       isSpatiallyReferenced = false;
       if isreal(cdata)
           hIm = images.internal.basicImageDisplay(hFig,hAx,...
                 squeeze(cdata(:,:,hImagePage.pageindex, :)),cdatamapping,clim,map,xdata,ydata,'nearest',isSpatiallyReferenced);
       else
           data_tmp = complex2color(cdata(:,:,hImagePage.pageindex), clim);
%            cdatamapping = 'direct';
           hIm = images.internal.basicImageDisplay(hFig,hAx,...
                 data_tmp,cdatamapping,clim,map,xdata,ydata,'nearest',isSpatiallyReferenced);
       end
       warning(warnstate);
       
       % listen page up/down
       pageListener = addlistener(hImagePage,'pageChanged',@imagePagedCData);
       setappdata(hIm,'updateListeners',pageListener);
       clear pageListener;
       
       % If it is an RSet base image tag it appropriately
       if ~isempty(rset)
           set(hIm,'tag','rset overview');
       end
       
       % Explicitly create an image model for the image.
       imgModel = getimagemodel(hIm);
       
       % Set original class type of imgmodel before image object is created.
       imgModel = setImageOrigClassType(imgModel,class(cdata));
       
       % For docked figures, we suppress the stream of warnings caused by
       % the figure positioning code in createPanels.
       if strcmpi(get(hFig,'WindowStyle'),'docked')
           old_state = warning('off','MATLAB:Figure:SetPosition');
           restoreWarningState = onCleanup(@() warning(old_state));
       end
       createPanels;
       clear restoreWarningState;
       
       % Enable all menus and toolbar buttons on list
       set(emptyToolInactiveList,'Enable','on');
       
       % attach button pressed callback for zoom in/out
       iptaddcallback(hFig,'WindowKeyPressFcn',@figureKeyPressed);
       
       % Enable Magnification combo box and wire to scrollpanel
       cb_api.setScrollpanel(sp_api);
       
       is_imtool_empty = false;
       
       if (~isSupportedImcontrastImage)
         set(contrastItems, 'Enable', 'off')
         
       else
         % Must be called after hIm defined
         imtoolModeManager.addMode(winLevelTool,...
                                   winLevelItem,...
                                   @makeWindowLevelModeCurrent,...
                                   @reactToModeChange);
       end
       
       if (~isSupportedImcolormap)
         set(colormapItem,'Enable','off');
       end
       
       % force a call to the resize function
       resizeImtool();
       
       set(hFig,'Visible','on'); 
       set(hFig,'ResizeFcn',@resizeImtool);
       
       % refresh image model if the image changes
       if ~isempty(hIm) && ishghandle(hIm)
           reactToImageChangesInFig(hIm,hFig,[],@refreshImageModel);
       end
       
     catch ME
         if ~strcmp(ME.identifier, 'images:getImageFromFile:OutOfMemTif')
             % In case of an unknown exception, close the figure.
             close(hFig)
         end
         rethrow(ME)

     end
     
     % Turn on HitTest so ButtonDownFcn will fire when image is clicked
     set(hIm,'HitTest','on')
     
     if ~isempty(rset)

         % Wire up image tiling.
         tile_manager = iptui.RSetTileManager(rset,hFig,sp_api);
         sp_api.addNewLocationCallback(@(pos) tile_manager.updateView());
         
         % Bug. In tile_manager.updateView we now need a way to force
         % updateView to happen even if button isn't down. Idea: pass an
         % additional argument to updateViewNewest(TF). updateView will
         % pass obj.buttonUp. Clients can call updateView(true) to force
         % update.
         tile_manager.updateView(true);

         % We need to disable imtool menus and toolbar buttons for modes
         % that do not support RSets yet.
         disableToolsForLargeImage(hFig);
         
     end
         
       %------------------------------------------------
       function common_args = parseFilename(common_args)
           
           if isrset(common_args.Filename)
               [common_args.CData,common_args.XData,common_args.YData,common_args.Map] = ...
                   getOverviewFromRSet(common_args.Filename);
           else
               
               try
                   [common_args.CData,common_args.Map] = ...
                       images.internal.getImageFromFile(common_args.Filename);
               catch ME
                   
                   out_of_mem_large_tif = strcmp(ME.identifier,'images:getImageFromFile:OutOfMemTif');
                   if out_of_mem_large_tif 
                     
                       rset_button = questdlg(getString(message('images:imtool:promptRsetCreation')),...
                                              getString(message('images:imtool:imageTooLarge')),...
                                              getString(message('images:commonUIString:create')),...
                                              getString(message('images:commonUIString:browse')),...
                                              getString(message('images:commonUIString:cancel')),...
                                              getString(message('images:commonUIString:create')));
                       
                       if isempty(rset_button) || ...
                          strcmpi(rset_button, getString(message('images:commonUIString:cancel')))

                           % User cancelled or hit "X" button when asked to
                           % create RSet.
                           rethrow(ME)
                           
                       elseif strcmpi(rset_button, getString(message('images:commonUIString:browse')))
                           
                           [common_args.CData,common_args.XData,common_args.YData,common_args.Map] = browseButtonPressed(ME);
                           
                       else
                         
                           [common_args.CData,common_args.XData,common_args.YData,common_args.Map] = yesButtonPressed(ME);
                           
                       end
                   else
                       rethrow(ME);
                   end
               end
           end
           
           %------------------------------------------------------
           function [cdata,xdata,ydata,map] = yesButtonPressed(ME)
                              
               repeat_dialog_display = true;
               while repeat_dialog_display
                   % User hit Yes button.
                   [temp,f_name] = fileparts(common_args.Filename); %#ok<ASGLU>
                   [f_name,pathname] = uiputfile({'*.rset','R-Set File'},...
                       'Save R-Set',...
                       fullfile(pwd, sprintf('%s.rset', f_name)));
                   
                   user_cancelled = f_name == 0;
                   if user_cancelled
                       rethrow(ME)
                   end
                   
                   full_name = fullfile(pathname,f_name);
                   
                   try
                       %Createrset will fail if user specifies
                       %non-writable directory.
                       rset_name = rsetwrite(common_args.Filename,full_name);
                       repeat_dialog_display = false;
                   catch rsetwriteException %#ok<NASGU>
                       % rsetwrite will fail if user chooses
                       % non-writable directory.
                       h_dlg = errordlg(getString(message('images:imtool:cannotSaveFile')),...
                           getString(message('images:imtool:cantWriteRset')),...
                           'modal');
                       uiwait(h_dlg);
                   end
                   
               end
               
               if ~isempty(rset_name)
                   % We succeeded at creating an RSet for a
                   % TIF that was too large to fit in memory.
                   [cdata,xdata,ydata,map] = ...
                       getOverviewFromRSet(rset_name);
               else
                   rethrow(ME)
               end
               
           end %yesButtonPressed
           
           %---------------------------------------------------------
           function [cdata,xdata,ydata,map] = browseButtonPressed(ME)
               
               repeat_dialog_display = true;
               while repeat_dialog_display
                   [f_name,pathname] = uigetfile({'*.rset','R-Set File';...
                                                  '*','All Files'},...
                                                  getString(message('images:imtool:openRset')));
                   
                   rset_name = fullfile(pathname,f_name);
                   
                   try
                       is_valid_rset = isrset(rset_name);
                   catch %#ok<CTCH>
                       % We reach this point if the
                       % user cancels the uigetfile dialog and
                       % failed to specify an rset. Rethrow
                       % exception that got us here.
                       rethrow(ME)
                   end
                   
                   if is_valid_rset
                       
                       [cdata,xdata,ydata,map] = ...
                           getOverviewFromRSet(rset_name);
                       repeat_dialog_display = false;
                   else
                       % We get here if user successfully selected
                       % a .rset file from uigetfile dialog that
                       % isn't a valid .rset file.
                       h_dlg = errordlg(getString(message('images:imtoolUIString:badRSetErrorDlgMessage')),...
                           getString(message('images:imtoolUIString:badRSetErrorDlgTitle')),'modal');
                       uiwait(h_dlg);                       
                   end
                   
               end
               
           end %browseButtonPressed
           
       end %parseFilename
            
       %---------------------------------
       function figureKeyPressed(obj,evt) %#ok<INUSL>

           % allow control +/- to zoom in/out of imtool
           control_key_pressed = isequal(numel(evt.Modifier),1) && ...
               strcmpi(evt.Modifier{1},'control');

           if control_key_pressed
               switch (evt.Key)
                   case 'add'
                       currentMag = sp_api.getMagnification();
                       newMag = images.internal.findZoomMag('in',currentMag);
                       sp_api.setMagnification(newMag)
                   case 'subtract'
                       currentMag = sp_api.getMagnification();
                       newMag = images.internal.findZoomMag('out',currentMag);
                       sp_api.setMagnification(newMag)
               end
           end
           
       end % figureKeyPressed
       
       %----------------------------------
       function refreshImageModel(obj,evt) %#ok<INUSD>
           
           % delete old model
           rmappdata(hIm,'imagemodel');
           
           % delete knowledge of origin filename if one existed
           filename = [];
           
           % reset figure name
%            imtool_name = getString(message('images:imtoolUIString:toolNameWithNumber',tool_number));
           imtool_name = sprintf('CT Tool %d - %s', tool_number,image_name);
           set(hFig,'Name',imtool_name);
           
           % Explicitly create an image model for the image.
           imgModel = getimagemodel(hIm);
           
           % Set original class type of imgmodel before image object is created.
           imgModel = setImageOrigClassType(imgModel,class(cdata));
           
       end % refreshImageModel

   end % addImageToImtool

    %---------------------------------------
    function switchModes = reactToModeChange
        
        %TODO: Explain what is going on in this function.
        switchModes = true;
        
        cropRect = getappdata(hFig,'imcropRectButtonDownOld');
        cropRectOverImage = ~isempty(cropRect) && isvalid(cropRect);
        
        if cropRectOverImage
            
            yesString = getString(message('images:commonUIString:yes'));
            noString = getString(message('images:commonUIString:no'));
            cancelString = getString(message('images:commonUIString:cancel'));
 
            buttonName = questdlg(getString(message('images:imtoolUIString:cropDlgMessage')),...
                getString(message('images:imtoolUIString:cropDlgTitle')),...
                yesString,...
                noString,...
                cancelString,...
                yesString);
            
            switch buttonName
                case (yesString)
                    cropRect.completeCrop();
                case (noString)
                    cropRect.delete();
                case (cancelString)
                    switchModes = false;
            end
            
        end
        
    end % reactToModeChange

   %-------------------------------------------
   function createPanels
     % This function adds the necessary panels. e.g.
     %  imdisplayrange, imscrollpanel, impixelinfo

    bottomPanelHeight = 21; % initialize to reasonable height
    
     hScrollPanel = imscrollpanel(hFig,hIm);
     sp_api = iptgetapi(hScrollPanel);
     
     if is_imtool_empty
         
       % If the image is being imported, we will use whichever is 
       % smaller: 100%  or fitMag      
       initial_mag = 100;
       if (sp_api.findFitMag() < 1)
         initial_mag = 'fit';
       end
       
     end
     
     hFigPos =  matlab.ui.internal.PositionUtils.getDevicePixelPosition(hFig);
     
     % Create bottom panel
     hBottomPanel = uipanel('Parent',hFig,...
                            'Tag','bottom panel',...
                            'BorderType','none',...
                            'Units','Pixels');
     matlab.ui.internal.PositionUtils.setDevicePixelPosition(hBottomPanel,...
         [1 1 hFigPos(3) bottomPanelHeight]);
    
     iptui.internal.setChildColorToMatchParent(hBottomPanel,hFig);
     
     % Create tools in bottom panel
     hPixInfoPanel = impixelinfo(hBottomPanel,hIm);
     
%      showDisplayRange = strcmp(getImageType(imgModel),'intensity');
     showDisplayRange = true;
     
     if showDisplayRange
       hRangePanel = myimdisplayrange(hBottomPanel,hIm,hImagePage);
     end
%      hRangePanel = myimdisplayrange(hBottomPanel,hIm,hImagePage);
     
     % Things may not settle down since the creation of hBottomPanel with
     % units in pixels and creation of hPixInfoPanel parented to
     % hBottomPanel. Hence the drawnow before querying the position in the
     % next line. (See g983968).
     drawnow;
     
     % Adjust height of bottom panel
     hPixInfoPos = matlab.ui.internal.PositionUtils.getDevicePixelPosition(hPixInfoPanel);
     bottomPanelHeight = hPixInfoPos(4);
     pos = matlab.ui.internal.PositionUtils.getDevicePixelPosition(hBottomPanel);
     pos(4) = bottomPanelHeight;
     matlab.ui.internal.PositionUtils.setDevicePixelPosition(hBottomPanel,pos);
     
     % Use getImWidth/getImHeight to get overall spatial extent of
     % image or rset.
     im_width  = images.internal.getImWidth(hIm);
     im_height = images.internal.getImHeight(hIm);
     
     % set up magnification and/or figure size based on input parm
     switch initial_mag
      case 'fit'
        % for the 'fit' case, the image has to fit within the default
        % figure, while retaining at least one pixel in each dimension.
        % Since the figure will not change size, take the bottom panel into
        % account.
         if isImageTooThin( im_width, im_height,  ...
                            hFigPos(3), (hFigPos(4) - bottomPanelHeight))

            issueExtremeAspectRatioWarning;
            mag = 1;
        else
            resizeImtool % force layout update, must be before call to findFitMag
            mag = sp_api.findFitMag();
        end
         sp_api.setMagnification(mag);
       
      case 'adaptive'
       % Try for 100%
       mag = 1;
       sp_api.setMagnification(mag)
       
       is_image_small = im_width<=minFigWidth && im_height<=minFigHeight;
       
       if is_image_small
         % Set figure size to minimum, and don't bother calling initSize
         pos = matlab.ui.internal.PositionUtils.getDevicePixelPosition(hFig);
         pos(3) = minFigWidth;
         pos(4) = minFigHeight;
         matlab.ui.internal.PositionUtils.setDevicePixelPosition(hFig,pos)
         
       else
         % Adaptively adjust figure size
         isBorderTight = true;
         wid = 'images:initSize:adjustingMag';
         state = warning('query',wid);
         warning('off',wid)
         images.internal.initSize(hIm,mag,isBorderTight)
         warning(state)

         hFigPos = matlab.ui.internal.PositionUtils.getDevicePixelPosition(hFig); %refresh height after initSize
          
         % find the magnification that fits in the new scroll panel size
         adaptiveMag = sp_api.findMagnification(im_width, im_height);

         % Save the height set by initSize for later use,
         % in determining  need to expand figure. That figure size is used
         % even if we decide to display at 100%. The value of 'adaptiveMag'
         % can be fractional. We want to define the minimum number of 
         % pixels required to represent the magnification of the image as  
         % an integral value. 
         onScreenImH = ceil(adaptiveMag * im_height);
         
         % check that the image can be shrunk to fit the figure size set by
         % initSize while maintaining at least one pixel in each dimension.  
         % The figure will expand later if needed, for the bottom panel, so
         % its height is not subtracted here.
         if isImageTooThin( im_width, im_height,  ...
                            hFigPos(3),  hFigPos(4))
             
             issueExtremeAspectRatioWarning;
             mag = 1;
         else
             mag = adaptiveMag;
         end
         sp_api.setMagnification(mag) ;
         
         % Figure out if figure needs to grow to fit bottom panel
         % Must be called after call to initSize and magnification               
         % decisions are completed
         heightNeeded = onScreenImH + bottomPanelHeight;
         needMoreRoomForBottomPanel = heightNeeded > hFigPos(4);
         
         if needMoreRoomForBottomPanel
           hFigPos(2) = hFigPos(2) - (heightNeeded - hFigPos(4));
           hFigPos(4) = heightNeeded;
           matlab.ui.internal.PositionUtils.setDevicePixelPosition(hFig,hFigPos);
           movegui(hFig,'onscreen');
         end
         
       end % end if is_image_small
       
       resizeImtool % force layout update for all adaptive cases           
       
     otherwise
   
         sp_api.setMagnification(initial_mag/100);
         
     end   % end of switch statement
     
        %-------------------------------------------
        function image_is_too_thin = isImageTooThin(...
                                                    imageWidth,...
                                                    imageHeight,...
                                                    destinationW,...
                                                    destinationH)
        %   This function compares the image to the rectangle 
        %   there it will be displayed.  If the aspect ratio is
        %   such that the image cannot be scaled to fit with each of its 
        %   dimensions being at least 1 pixel, this function returns true.

           image_is_too_thin =  imageWidth/imageHeight > destinationW ||...
                                imageHeight/imageWidth > destinationH;
       end
       
       %--------------------------------------
       function issueExtremeAspectRatioWarning
       %    The function issues a warning that the image cannot be
       %    displayed as requested with at least one pixel showing in
       %    each dimension
          
           warning(message('images:imtool:extremeAspectRatio'));
           
       end   
       
   end % end createPanels
     
   %---------------------------------------
   function toolbar = createMenusAndToolbar

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % File menu 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     fileMenu = uimenu(hFig,...
         'Label',getString(message('images:imtoolUIString:fileMenubarLabel')),...
         'Tag','file menu');

     % Common properties
     m.Parent = fileMenu;
     
     % New item
     m.Label       = getString(message('images:imtoolUIString:newMenubarLabel'));
     m.Accelerator = 'N';
     m.Callback    = @openNewTool;
     m.Tag         = 'new menu item';
     m.Separator   = 'off';     
     uimenu(m)
   
     % Open item
     m.Label       = getString(message('images:imtoolUIString:openMenubarLabel'));
     m.Accelerator = 'O';
     m.Callback    = @openFromFile;
     m.Tag         = 'open menu item';
     m.Separator   = 'off';     
     uimenu(m)
     
     % Import from Workspace item
     m.Label       = getString(message('images:imtoolUIString:importMenubarLabel'));
     m.Accelerator = '';
     m.Callback    = @openFromWS;
     m.Tag         = 'open from workspace menu item';
     m.Separator   = 'off';     
     uimenu(m)

     % Export to Workspace item
     m.Label       = getString(message('images:imtoolUIString:exportMenubarLabel'));
     m.Accelerator = 'E';
     m.Callback    = @callImExportToWorkspace;
     m.Tag         = 'export to workspace menu item';     
     m.Separator   = 'on';     
     exportToWSItem = uimenu(m);
     
     % Save As item
     m.Label       = getString(message('images:imtoolUIString:saveAsMenubarLabel'));
     m.Accelerator = 'S';
     m.Callback    = @saveImage;
     m.Tag         = 'save as menu item';     
     m.Separator   = 'off';     
     saveAsItem = uimenu(m);
     
     % Preferences item
     m.Label       = getString(message('images:imtoolUIString:preferencesMenubarLabel'));
     m.Accelerator = '';
     m.Callback    = @(varargin) preferences('Image Processing');
     m.Tag         = 'preferences menu item';
     m.Separator   = 'on';
     uimenu(m);
     
     % Print to Figure item
     m.Label       = getString(message('images:imtoolUIString:printToFigureMenubarLabel'));
     m.Accelerator = '';
     m.Callback    = @printToFig;
     m.Tag         = 'print to figure menu item';
     m.Separator   = 'on';
     printToFigItem = uimenu(m);
     
     % Close item
     m.Label       = getString(message('images:imtoolUIString:closeMenubarLabel'));
     m.Accelerator = 'W';
     m.Callback    = @closeImtool;
     m.Tag         = 'close menu item';
     m.Separator   = 'on';
     uimenu(m)
     
     % Close all item
     m.Label       = getString(message('images:imtoolUIString:closeAllMenubarLabel'));
     m.Accelerator = '';
     m.Callback    = @closeAll;
     m.Tag         = 'close all menu item';
     m.Separator   = 'off';
     uimenu(m)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % end of File menu creation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % Tools menu 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     toolsMenu = uimenu(hFig,...
         'Label',getString(message('images:imtoolUIString:toolsMenubarLabel')),...
         'Tag','tools menu');
     
     % Common properties for tools menu items
     m.Parent      = toolsMenu;

     % Add mode menu items first because they are used more than
     % action menu items.
     % Note: This is different from the toolbar which puts mode items last.
     
     % zoom in menu item
     m.Label       = getString(message('images:commonUIString:zoomInMenubarLabel'));
     m.Tag         = 'zoom in menu item';
     zoomInItem = uimenu(m);
 
     % zoom out menu item
     m.Label       = getString(message('images:commonUIString:zoomOutMenubarLabel'));
     m.Tag         = 'zoom out menu item';
     zoomOutItem = uimenu(m);
 
     % pan menu item
     m.Label       = getString(message('images:commonUIString:panMenubarLabel'));
     m.Tag         = 'pan menu item';
     panItem = uimenu(m);
      
     % windowlevel menu item
     m.Label       = getString(message('images:imtoolUIString:windowLevelMenubarLabel'));
     m.Tag         = 'windowlevel menu item';
     winLevelItem = uimenu(m);
 
     % crop tool menu item
     m.Label       = getString(message('images:imtoolUIString:cropImageMenubarLabel'));
     m.Separator   = 'on';
     m.Tag         = 'crop tool menu item';
     cropItem = uimenu(m);
     
     % distance tool menu item
     m.Label       = getString(message('images:imtoolUIString:measureDistanceMenubarLabel'));
     m.Separator   = 'off';
     m.Tag         = 'distance tool menu item';
     distanceItem = uimenu(m);
     
     % overview menu item
     m.Label       = getString(message('images:imtoolUIString:overviewMenubarLabel'));
     m.Callback    = @showOverview;
     m.Tag         = 'overview menu item';
     m.Separator   = 'on';     
     overviewItem = uimenu(m);
     
     % pixel region menu item
     m.Label       = getString(message('images:imtoolUIString:pixelRegionMenubarLabel'));
     m.Separator   = 'off';
     m.Callback    = @showPixelRegionTool;
     m.Tag         = 'pixel region menu item';
     pixelRegionItem = uimenu(m);
     
     % image info menu item
     m.Label       = getString(message('images:imtoolUIString:imageInformationMenubarLabel'));
     m.Separator   = 'off';
     m.Callback    = @showImageInfo;
     m.Tag         = 'image info menu item';
     infoItem = uimenu(m);

     % adjust contrast menu item
     m.Label       = getString(message('images:imtoolUIString:adjustContrastMenubarLabel'));
     m.Separator   = 'off';
     m.Callback    = @showImcontrast;
     m.Tag         = 'adjust contrast menu item';
     adjustContrastItem = uimenu(m);
     
     % choose colormap menu button
     m.Label       = getString(message('images:imtoolUIString:chooseColormapMenubarLabel'));
     m.Separator   = 'on';
     m.Callback    = @showImColormapSelect;
     m.Tag         = 'choose colormap menu item';
     colormapItem = uimenu(m);
     
%      % page down menu button
%      m.Label       = 'page down';
%      m.Separator   = 'on';
%      m.Callback    = @(x) []; % do nothing
%      m.Tag         = 'page down (in working)';
%      pagedownItem  = uimenu(m);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % end of Tools menu creation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % Window menu
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
     matlab.ui.internal.createWinMenu(hFig);


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % end of Window menu creation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % Help menu
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
     helpMenu = uimenu(hFig,...
         'Label',getString(message('images:imtoolUIString:helpMenubarLabel')),...
         'Tag','help menu');

     % imtool help menu item
     m.Parent   = helpMenu;
     m.Label    = getString(message('images:imtoolUIString:imageToolHelpMenubarLabel'));
     m.Callback = @showHelp;
     m.Tag      = 'help menu item';
     uimenu(m)
     
     % ipt standard help menu items
     iptstandardhelp(helpMenu);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % end of Help menu creation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Toolbar
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     toolbar =  uitoolbar(hFig);
 
     % Get roots for where to find icons
     [iconRoot, iconRootMATLAB] = ipticondir;
     iconRootMytools = fullfile(fileparts(which('CTsimulation')), 'tools\cttool\images');
     
     % Common properties for several push tools
     t = [];
     t.toolConstructor            = @uipushtool;
     t.properties.Parent          = toolbar;
     t.properties.Interruptible   = 'off';    % Make buttons busy, cancel
     t.properties.BusyAction      = 'cancel'; % any repeated click events.
      
     % overview toolbar button
     t.iconConstructor            = @makeToolbarIconFromPNG;    
     t.iconRoot                   = iconRoot;    
     t.icon                       = 'overview.png';
     t.properties.ClickedCallback = @showOverview;
     t.properties.TooltipString   = getString(message('images:imtoolUIString:overviewTooltipString'));
     t.properties.Tag             = 'overview toolbar button';
     overviewTool = makeToolbarItem(t);
      
     % pixel region toolbar button
     t.iconConstructor            = @makeToolbarIconFromPNG;    
     t.iconRoot                   = iconRoot;    
     t.icon                       = 'pixel_region.png';
     t.properties.ClickedCallback = @showPixelRegionTool;
     t.properties.TooltipString   = getString(message('images:imtoolUIString:pixelRegionTooltipString'));
     t.properties.Tag             = 'pixel region toolbar button';
     pixelRegionTool = makeToolbarItem(t);
         
     % image info toolbar button
     t.iconConstructor            = @makeToolbarIconFromPNG;      
     t.iconRoot                   = iconRoot;    
     t.icon                       = 'icon_info.png';
     t.properties.ClickedCallback = @showImageInfo;
     t.properties.TooltipString   = getString(message('images:imtoolUIString:imageInfoTooltipString'));
     t.properties.Tag             = 'image info toolbar button';
     infoTool = makeToolbarItem(t);
     
     % adjust contrast toolbar button
     t.iconConstructor            = @makeToolbarIconFromPNG;
     t.iconRoot                   = iconRoot;    
     t.icon                       = 'tool_contrast.png';
     t.properties.ClickedCallback = @showImcontrast;
     t.properties.TooltipString   = getString(message('images:imtoolUIString:adjustContrastTooltipString'));
     t.properties.Tag             = 'adjust contrast toolbar button';
     adjustContrastTool = makeToolbarItem(t);
      
     % help toolbar button
     t.iconConstructor            = @makeToolbarIconFromGIF;
     t.iconRoot                   = iconRootMATLAB;    
     t.icon                       = 'helpicon.gif';
     t.properties.ClickedCallback = @showHelp;
     t.properties.TooltipString   = getString(message('images:commonUIString:help'));
     t.properties.Tag             = 'imtool help toolbar button';
     makeToolbarItem(t);
          
     % crop tool toolbar button
     t.toolConstructor            = @uitoggletool;
     t.properties.Parent          = toolbar;
     t.iconConstructor            = @makeToolbarIconFromPNG;
     t.iconRoot                   = iconRoot;
     t.icon                       = 'crop_tool.png';
     t.properties.TooltipString   = getString(message('images:imtoolUIString:cropImageTooltipString'));
     t.properties.Tag             = 'crop tool toolbar button';
     cropTool = makeToolbarItem(t);
     
     % distance tool toolbar button
     t.toolConstructor            = @uitoggletool;
     t.properties.Parent          = toolbar;
     t.iconConstructor            = @makeToolbarIconFromGIF;
     t.iconRoot                   = iconRoot;
     t.icon                       = 'distance_tool.gif';
     t.properties.TooltipString   = getString(message('images:imtoolUIString:measureDistanceTooltipString'));
     t.properties.Tag             = 'distance tool toolbar button';
     distanceTool = makeToolbarItem(t);
     
     % adjust pageup toolbar button
     t.toolConstructor            = @uipushtool;
     t.iconConstructor            = @makeToolbarIconFromPNG;
     t.iconRoot                   = iconRootMytools;    
     t.icon                       = 'leftarrow.png';
     t.properties.ClickedCallback = @imagePageup;
     t.properties.TooltipString   = 'page up';
     t.properties.Tag             = 'page up toolbar button';
     PageupTool = makeToolbarItem(t);
     
     % adjust pagedown toolbar button
     t.toolConstructor            = @uipushtool;
     t.iconConstructor            = @makeToolbarIconFromPNG;
     t.iconRoot                   = iconRootMytools;    
     t.icon                       = 'rightarrow.png';
     t.properties.ClickedCallback = @imagePagedown;
     t.properties.TooltipString   = 'page down';
     t.properties.Tag             = 'page down toolbar button';
     PagedownTool = makeToolbarItem(t);
     
     
     
     % Create mode toolbar items at end of toolbar
     
     % navigational toolbar buttons (zoom in, zoom out, pan)
     navToolButtons = navToolFactory(toolbar);
     
     % Put each tool button in its variable that has function scope so modes can be
     % set up correctly in the main function.
     zoomInTool = navToolButtons.zoomInTool;
     zoomOutTool = navToolButtons.zoomOutTool;
     panTool = navToolButtons.panTool;
     
     % Place separator before zoom in toolbar button
     set(navToolButtons.zoomInTool,'Separator','on')
     
     % Place separator before crop tool toolbar button
     set(cropTool,'Separator','on');
     
     % window level toolbar button
     t = []; % reset structure 
     t.toolConstructor            = @uitoggletool;
     t.iconConstructor            = @makeToolbarIconFromPNG;
     t.iconRoot                   = iconRoot;    
     t.icon                       = 'cursor_contrast.png';
     t.properties.Parent          = toolbar;
     t.properties.TooltipString   = getString(message('images:imtoolUIString:windowLevelTooltipString'));
     t.properties.Tag             = 'windowlevel toolbar button';
     winLevelTool = makeToolbarItem(t);
     

     
     % We must ensure that all pending layout of toolbar is complete before
     % magComboBox layout begins, otherwise layout of magComboBox can
     % happen before other toolbar elements are added. g929824
     drawnow('expose')
    
     % magnification combo box
     [cb, cb_api] = createMagComboBox(toolbar); %#ok<ASGLU>
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % end of toolbar creation
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
     % Make lists of items that need to be enabled or disabled as a group
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
     
     % This is a list of handles to toolbar and menu items to disable when the adjust
     % contrast and window level tools will not work on an image.
     contrastItems = [winLevelTool, winLevelItem,...
                      adjustContrastTool, adjustContrastItem];
 
     % This is a list of handles to toolbar buttons and menu items to disable when
     % imtool is opened in an empty state.
     emptyToolInactiveList = [
         zoomInTool, zoomOutTool, panTool, winLevelTool,...
         overviewTool, infoTool, distanceTool, cropTool,...
         pixelRegionTool, adjustContrastTool,PageupTool,PagedownTool, ...
         zoomInItem, zoomOutItem, panItem, winLevelItem,...
         printToFigItem, colormapItem, ...
         exportToWSItem, saveAsItem, overviewItem, pixelRegionItem,...
         distanceItem, cropItem, infoItem,adjustContrastItem];
              
   end
      
   %-----------------------------------------------  
   function callImExportToWorkspace(varargin)
     
     imExportToWorkspace(hFig);
     
   end %callImExportToWorkspace

   %-----------------------------------------------  
   function saveImage(varargin)
     imsave(hFig);
   end %saveImage
   
  %-----------------------------------------------
  function closeImtool(varargin)
    close(hFig);
  end
  
  %-----------------------------------------------
  function openNewTool(varargin)
    imtool;
  end

  %-----------------------------------------------
  function openFromFile(varargin)
    [filename,user_canceled] = imgetfile;
    if user_canceled
      return
    else
      if is_imtool_empty
          try
            addImageToImtool(filename,'initialmag','fit');
          catch ME
              % If addImageToImtool fails with a TIF out of memory
              % exception, swallow error. Otherwise, rethrow exception.
              if ~strcmp(ME.identifier,'images:getImageFromFile:OutOfMemTif')
                  rethrow(ME)
              end
          end  
      else
        imtool(filename);
      end
    end
  end

  %-----------------------------------------------
  function openFromWS(varargin)
    [I,map,I_var_name,map_var_name,user_canceled] = iptui.internal.myimgetvar(hFig);
      
    if user_canceled
       return
    end
    
    if is_imtool_empty            
      image_name = I_var_name;
      if isempty(map)
        addImageToImtool(I,'initialmag','fit');
      else
        addImageToImtool(I,map,'initialmag','fit');
      end
    else
        isImageTypeIndexed = ~isempty(map);
        if isImageTypeIndexed
            if isempty(map_var_name)
                % user is importing an indexed image, but did not select a
                % colormap.
                evalin('base',sprintf('imtool(%s,gray(256))',I_var_name));
            else
                evalin('base',sprintf('imtool(%s,%s)',I_var_name,map_var_name));
            end
        else
            evalin('base',sprintf('imtool(%s);',I_var_name));
        end
    end
    
  end

  %-----------------------------------------------
  function printToFig(varargin)
    
    printImageToFigure(hScrollPanel);
    
  end
  
  %-----------------------------------------------
  function closeAll(varargin)
     figs = findall(0,'Type','figure','Tag','imtool');
     close(figs)
  end
    
  %------------------------------
  function showOverview(varargin)
  
    hOverviewFig = showChildTool(hFig,hOverviewFig, @imoverview, {hIm},...
                                 getString(message('images:commonUIString:overview')),...
                                 tool_number,...
                                 {'left','right'},...
                                 {'top','top'});
  end

  %-------------------------------
  function showImageInfo(varargin)

    if isempty(filename)
        args = {hIm};
    elseif ~isempty(rset)
        args = {hIm, rset.getRSetDetails()};
    else
        args = {hIm, filename};
    end
  
    hImageInfoFig = showChildTool(hFig,hImageInfoFig, @imageinfo, args,...
                                 getString(message('images:commonUIString:imageInformation')),...
                                  tool_number,...
                                  {'left','left'},...
                                  {'bottom','top'});
  end

  %--------------------------------
  function showImcontrast(varargin)
      
      if (isempty(rset))
          hContrastFig = showChildTool(hFig,hContrastFig, @myimcontrast, {hIm},...
                                         getString(message('images:commonUIString:adjustContrast')),...
                                         tool_number,...
                                             {'left','left'},...
                                         {'bottom','top'});
      else
          rsetMode = true;
          hContrastFig = showChildTool(hFig,hContrastFig, @myimcontrast, {hIm, rsetMode},...
                                         getString(message('images:commonUIString:adjustContrast')),...
                                         tool_number,...
                                             {'left','left'},...
                                         {'bottom','top'});
      end
      
      if isempty(hContrastFig) || ~ishghandle(hContrastFig)
          return
      end

  end
  
  %-------------------------------------
  function showPixelRegionTool(varargin)

    hPixelRegionFig = showChildTool(hFig,hPixelRegionFig, ...
                                    @impixelregion, {hIm},...
                                    getString(message('images:commonUIString:pixelRegion')),...
                                    tool_number,...
                                    {'right','left'},...
                                    {'top','top'});

    % align pixelregion view with imtool view
    pixregion_sp = findobj(hPixelRegionFig,'Tag','imscrollpanel');
    pixregion_api = getappdata(pixregion_sp,'impixelregionpanelAPI');
    pixregion_api.centerRectInViewport();
    
  end
  
  %-------------------------------------
  function showImColormapSelect(varargin)
    
    hColormapSelectFig = showChildTool(hFig,hColormapSelectFig, ...
                                       @imcolormaptool, {hFig},...
                                       getString(message('images:commonUIString:chooseColormap')),...
                                       tool_number,...
                                       {'right','left'},...
                                       {'top','top'});
  end
    
  %-----------------------------------
  function showHelp(varargin)
  
    topic = 'imageframe';
    helpview([docroot '/toolbox/images/images.map'],topic);

  end    
  
  %------------------------------------
  function imagePagedown(varargin)
    hImagePage.imagePageDown();
  end

  %------------------------------------
  function imagePageup(varargin)
    hImagePage.imagePageUp();
  end
  
  %------------------------------------
  function imageScroll(varargin)
    eventdata = varargin{2};
    if eventdata.VerticalScrollCount > 0
        % page down
        imagePagedown;
    elseif eventdata.VerticalScrollCount < 0
        % page up
        imagePageup;
    end
  end
  %------------------------------------
  function imagePagedCData(varargin)
      if isreal(cdata)
          set(hIm, 'CData', squeeze(cdata(:,:,hImagePage.pageindex, :)));
      else
          data_tmp = complex2color(cdata(:,:,hImagePage.pageindex), clim);
          set(hIm, 'CData', data_tmp);
      end
  end

  %--------------------------------------
  function setCurrentMode(fun,ptr)
  % ptr is a cell array containing relevant pointer param/value
  % pairs.  These are param/value pairs returned by the SETPTR fcn.

    removeCurrentMode()    

    enterFcn = @(f,cp) set(f, ptr{:});
    iptSetPointerBehavior(hIm, enterFcn);
    sp_api.setImageButtonDownFcn(fun)
    
  end

  %-------------------
  function removeCurrentMode
    
    sp_api.setImageButtonDownFcn([])
    iptSetPointerBehavior(hIm, []);

  end

  %----------------------------------------
  function makeDefaultModeCurrent(varargin)

    removeCurrentMode()

  end

  %------------------------------------
  function makeCropModeCurrent(varargin)

    fun = @iptui.imcropRectButtonDown;
    ptr = setptr('crosshair');
    
    setCurrentMode(fun,ptr)
    
  end

  %------------------------------------
  function makeDistanceModeCurrent(varargin)

    fun = @iptui.imdistlineButtonDown;
    ptr = setptr('crosshair');
    
    setCurrentMode(fun,ptr)

  end

  %---------------------------------------
  function makeZoomInModeCurrent(varargin)

    fun = @imzoomin;
    ptr = setptr('glassplus');
    
    setCurrentMode(fun,ptr)

  end

  %----------------------------------------
  function makeZoomOutModeCurrent(varargin)

    fun = @imzoomout;
    ptr = setptr('glassminus');
    
    setCurrentMode(fun,ptr)

  end

  %------------------------------------
  function makePanModeCurrent(varargin)

    fun = @impan;
    ptr = setptr('hand');
    
    setCurrentMode(fun,ptr)

  end

  %--------------------------------------------
  function makeWindowLevelModeCurrent(varargin)

    fun = @(obj,evt)(windowlevel(hIm, hFig));
    ptr = {'Pointer', 'custom',...
           'PointerShapeCData', getWLPointer,...
           'PointerShapeHotSpot', [7 7]};
      
    setCurrentMode(fun,ptr)

  end

  %------------------------------  
  function resizeImtool(varargin)
  
    figPos = matlab.ui.internal.PositionUtils.getDevicePixelPosition(hFig);
    is_window_docked = strcmpi('docked',get(hFig,'WindowStyle'));
    
    if ~is_window_docked
      % Constrain minimum figure size
      if ~isempty(lastFigPos)
        figPos = constrainFigSize(figPos,lastFigPos,[minFigWidth minFigHeight]);
      end
      
      lastFigPos = figPos;
      matlab.ui.internal.PositionUtils.setDevicePixelPosition(hFig,figPos);
    end
    
    if showDisplayRange
        % Make sure hRangePanel moves L/R as window shrinks/grows.
        rangePos = matlab.ui.internal.PositionUtils.getDevicePixelPosition(hRangePanel);
        rangePos(1) = figPos(3)-rangePos(3);
        matlab.ui.internal.PositionUtils.setDevicePixelPosition(hRangePanel,rangePos);
    end

    % Update hBottomPanel width to match figure width
    bottomPanelPos = matlab.ui.internal.PositionUtils.getDevicePixelPosition(hBottomPanel);
    bottomPanelPos(3) = figPos(3);
    matlab.ui.internal.PositionUtils.setDevicePixelPosition(hBottomPanel,bottomPanelPos);    

    % Set hScrollPanel position in pixels to allow room for bottom panel
    spB = bottomPanelPos(4)+1;
    spH = figPos(4) - bottomPanelPos(4); 
    if (spH > 0) && (figPos(3) > 0)
        % set throws an error if Width or Height is zero or less
        set(hScrollPanel,'Units','pixels');
        matlab.ui.internal.PositionUtils.setDevicePixelPosition(hScrollPanel,...
            [1 spB figPos(3) spH]);
    end
    
  end
  
  %------------------  
  function deleteTool(varargin)
  % Delete all child windows that the image tool may have created.
  
    id_stream.recycleId(tool_number);
    
    if ~isempty(hOverviewFig) && ishghandle(hOverviewFig)
        delete(hOverviewFig)
    end
    
    if ~isempty(hPixelRegionFig) && ishghandle(hPixelRegionFig)
        delete(hPixelRegionFig);
    end

    if ~isempty(hImageInfoFig) && ishghandle(hImageInfoFig)
        delete(hImageInfoFig);
    end

    if ~isempty(hContrastFig) && ishghandle(hContrastFig)
        delete(hContrastFig);
    end

  end

  %--------------------------------------------
  function tf = isSupportedImcontrastImage
  
    if (~isequal(getImageType(imgModel), 'intensity'))
        tf = false;
        return
    end
    
    switch (getClassType(imgModel))
      case {'uint8', 'uint16', 'uint32'}
        tf = true;
      
      case {'int8', 'int16', 'int32'}
        tf = true;
      
      case {'logical'}
        tf = false;
      
      case {'double', 'single'}
        tf = true;
      
      otherwise
        tf = false;
      
    end
  
  end
  
  %--------------------------------------------        
  function tf = isSupportedImcolormap
  
    tf = ~strcmpi(getImageType(imgModel),'truecolor');
    
  end

  %---------------------------------------------------------------
  function [cdata,xdata,ydata,map] = getOverviewFromRSet(filename)

      % Cache rset at function scope. TODO: HOW TO DO THIS WITHOUT TWO
      % WARNINGS if base TIF not located.
      rset = iptui.RSet(filename);
      
      [cdata,xdata,ydata] = rset.getOverview();
      map = rset.getColormap();
    
  end

end %imtool

%------------------------------------------------------------------
function hChildFig = showChildTool(hParentFig,hChildFig,...
                                   childFun,childArgs,...
                                   childToolName,...
                                   parentToolNumber,...
                                   parentChildLocHor,...
                                   parentChildLocVer)

  if isempty(hChildFig) || ~ishghandle(hChildFig)
      hChildFig = feval(childFun,childArgs{:});
      
      if isempty(hChildFig) || ~ishghandle(hChildFig)
          return
      end

      childName = sprintf('%s (%s)', ...
                          childToolName, getString(message('images:imtoolUIString:toolNameWithNumber',parentToolNumber)));
                      
      set(hChildFig,'Name',childName)
      
      iptwindowalign(hParentFig, parentChildLocHor{1}, ...
                    hChildFig, parentChildLocHor{2});
      iptwindowalign(hParentFig, parentChildLocVer{1},...
                    hChildFig, parentChildLocVer{2});
      
  else
      figure(hChildFig);
  end
    
end

%------------------------------------------------------------------
function PointerShapeCData = getWLPointer

iconRoot = ipticondir;
cdata = makeToolbarIconFromPNG(fullfile(iconRoot,'cursor_contrast.png'));
PointerShapeCData = cdata(:,:,1)+1;

end

%-----------------------------------------
function width  = getToolbarWidth(toolbar)
  
  % This function is designed to calculate the toolbar width based on the number
  % of buttons and separators.  It is assumed there is one magnification combo
  % box.

  % The button icons are 16x16 pixels.  The other object sizes were measured
  % empirically using a figure of known width as a ruler.
  buttonWidth = 16;
  buttonBuffer = 8;
  separatorWidth = 8;
  magpanelWidth = 140;
  largeFontsBuffer = 50;
  
  % The magnification combobox javacomponent is currently not a child of the
  % toolbar.  This may change in the future.
  numButtons = numel(get(toolbar,'Children'));
  numSeparators = numel(findall(toolbar,'Separator','on'));
    
  width = numButtons*buttonWidth + numButtons*buttonBuffer + ...
        numSeparators*separatorWidth + magpanelWidth + largeFontsBuffer;
end

function disableToolsForLargeImage(hFig)

menuOptionsToDisable = ...
    findobj(hFig,'tag','save as menu item','-or',...
    'tag','export to workspace menu item','-or',...
    'tag','crop tool menu item','-or',...
    'tag','pixel region menu item','-or',...
    'tag','choose colormap menu item');

toolbarButtonsToDisable = ...
    findobj(hFig,...
    'tag','pixel region toolbar button','-or',...
    'tag','crop tool toolbar button');

set(menuOptionsToDisable,'Enable','off');
set(toolbarButtonsToDisable','Enable','off');

end


function dataout = complex2color(datain, clim)
    dataout = reshape(repmat(datain, 1, 1, 3), [], 3);
    dataout(:, 3) = (real(datain(:)) - clim(1))./(clim(2)-clim(1));
    dataout(dataout(:, 3)<0, 3) = 0;
    dataout(dataout(:, 3)>1, 3) = 1;
    dataout(:,2) = imag(datain(:)).*(4/(clim(2)-clim(1)));
    dataout(:,1) = (dataout(:,2)<0).*0.5 + 0.04;
    dataout(:,2) = abs(dataout(:,2));
    dataout(dataout(:,2)>1, 2) = 1;

    dataout = hsv2rgb(dataout);
    dataout = reshape(dataout, [size(datain) 3]);

end