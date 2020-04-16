function [dataflow, prmflow, status] = reconnode_CRISAxialBackprojection(dataflow, prmflow, status)
% external recon node, Axial FBP call CRIS FBP 
% [dataflow, prmflow, status] = reconnode_CRISAxialFBP(dataflow, prmflow, status);
% hard code for temporay using

% parameters set in pipe
FBPprm = prmflow.pipe.(status.nodename);

% parameters for recon
Nshot = prmflow.recon.Nshot;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nviewprot = prmflow.recon.Nviewprot;
Nview = prmflow.recon.Nview;

if isfield(FBPprm, 'FOV')
    reconFOV = FBPprm.FOV;
else
    reconFOV = prmflow.external.rawxml.ReconParameters.displayFOV;
end
if isfield(FBPprm, 'imagesize')
    imagesize = BPprm.imagesize;
else
    imagesize = 512;
end
if isfield(FBPprm, 'Kernel')
    reconKernel = FBPprm.Kernel;
else
    reconKernel = prmflow.external.rawxml.ReconParameters.reconKernel;
end

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

% GPU engine
Engine = EngineManager.GetInstance('GPUCompute');

% ini image
dataflow.image = zeros(imagesize, imagesize, Nslice*Nshot, 'single');

% loop the shots
for ishot = 1:Nshot
    viewindex = (1:Nviewprot) + (ishot-1)*Nviewprot;
    struct_raw_Inside.data = dataflow.rawdata(:, :, viewindex);
    struct_raw_Inside.raw_size = [Npixel, Nslice, Nviewprot];
%     struct_raw_Inside.header.viewAngle = -90;   % 正反投影中所用到的库不是采用同一个坐标系，导致需要逆时针偏转90°

    % BP
    BpStruct.nRotateDirection    = -1;                                       % 投影扫描方向，1:clockwise, -1:counterclockwise
    BpStruct.nViewPerRevolution  = Nviewprot;                               % 投影采样个数
    BpStruct.nTotalViewNumber    = Nviewprot;                               % ?
    BpStruct.fStartAngle         = prmflow.recon.startviewangle(ishot) - pi/2;   
                                                                            % 投影起始角度
    BpStruct.nChannelDirection   = -1;                                      % 探测器排布方向，-1:clockwise, 1:counterclockwise
    BpStruct.nChannelNumPar      = Npixel;                                  % 探测器通道数
    BpStruct.fChannelParSpace    = prmflow.recon.delta_d;
    BpStruct.fMidChannelPar      = prmflow.recon.midchannel-1;                % 探测器中心通道 (-1)
    BpStruct.nSliceDirection     = 1;
    BpStruct.nSliceNumber        = Nslice;                
    BpStruct.fSliceSpace         = prmflow.system.detector.hz_ISO;
    % BpStruct.nCouchDirection     = gParas.AcquisitionParameter.tableDirection; % Add for 3DBP
    % BpStruct.fSourseToIsocenter  = gParas.GeometryParameter.sourceToIso; % Add for 3DBP
    BpStruct.fTiltAngle          = 0;
    BpStruct.fZDFS               = 0; % zDFS 情况下光源抖动偶数View位置相对奇数View位置的距离
    BpStruct.fMaxFOV             = 500;
    BpStruct.fReconFOV           = reconFOV;
    BpStruct.nXPixels            = imagesize;
    BpStruct.nYPixels            = imagesize;
    BpStruct.fXReconCenter       = 0;
    BpStruct.fYReconCenter       = 0;
    % BpStruct.fImageThickness     = gParas.ReconParameters.ImageThickness; % Add for 3DBP, NO Using
    % BpStruct.fImageIncrement     = gParas.ReconParameters.ImageIncrement; % Add for 3DBP
    BpStruct.nImageNumber        = BpStruct.nSliceNumber;

    fViewWeight = ones(BpStruct.nTotalViewNumber, 1, 'single')*0.5;
    ImgNum_pZ = BpStruct.nSliceNumber;
    
    func = Engine.GetFunction('Axial2DBP_CU');
    CorrImage = func(struct_raw_Inside.data, fViewWeight, BpStruct, ImgNum_pZ);
%     % rot back
%     CorrImage=rot90(CorrImage,1);
    % copy to dataflow, /2 to compatible with matlab iradon
    pageindex = (1:Nslice) + Nslice*(ishot-1);
    dataflow.image(:,:,pageindex) = permute(CorrImage, [2, 1, 3])./2;
end

%% done
% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end