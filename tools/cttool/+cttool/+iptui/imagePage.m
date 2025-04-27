classdef imagePage < handle
    properties
        Npage;
        pageindex;
    end
    events
        pageChanged;
    end
    methods (Access = public)
        % ini
        function obj = imagePage(Npage, pageindex)
            if nargin<2
                pageindex = 1;
            end
            obj.Npage = Npage;
            obj.pageindex = pageindex;
        end
        
        % page down
        function r = imagePageDown(obj)
            if obj.pageindex < obj.Npage
                obj.pageindex = obj.pageindex + 1;
                obj.notify('pageChanged');
                r = true;
            else
                r = false;
            end
        end
        
        % page uo
        function r = imagePageUp(obj)
            if obj.pageindex > 1
                obj.pageindex = obj.pageindex - 1;
                obj.notify('pageChanged');
                r = true;
            else
                r = false;
            end
        end
        
        % page reset
        function imagePageReset(obj, Npage, pageindex)
            if nargin<3
                pageindex = 1;
            end
            obj.Npage = Npage;
            obj.pageindex = pageindex;
            obj.notify('pageChanged');
        end
    end
end