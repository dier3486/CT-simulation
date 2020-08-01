

hFig = figure('Toolbar','none', ...
              'Tag','CTtool', ...
              'DeleteFcn',@deleteTool);

id_stream = idStreamFactory('ImtoolInstance');
tool_number = id_stream.nextId();
          
          
          
          
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
