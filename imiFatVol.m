classdef imiFatVol
    %IMIFATVOL Wrapper class for Imiomics fat volumes and data
    
    properties
        % numbers
        poemid % int           
        female % logical
        dxaval % float
        % 3D volumes
        DefVol
        JacDet
    end
    
    methods
        function obj = imiFatVol(arg1,arg2,arg3,arg4,arg5)
            %IMIFATVOL Construct an instance of this class
            % allow 0 args for initialization
            if nargin == 5
            obj.poemid = arg1;        
            obj.female = arg2;
            obj.dxaval = round(arg3);
            obj.DefVol = arg4;
            obj.JacDet = arg5;
            end
        end
        
        % constructor without volumes (set them later)
        function obj = imiFatVolNoVol(arg1,arg2,arg3)
            obj.poemid = arg1;   
            obj.female = arg2;
            obj.dxafat = arg3;
        end
        
         
    end
end

