function [dataflow, prmflow, status] = reconnode_CRISAxialFBP(dataflow, prmflow, status)
% external recon node, Axial FBP call CRIS FBP 
% [dataflow, prmflow, status] = reconnode_CRISAxialFBP(dataflow, prmflow, status);
% hard code for temporay using

% BP prepare
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

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

% Filter
% copy the rawxml to gParas_Inside
gParas_Inside = prmflow.external.rawxml;
% set the values,
gParas_Inside.VariableParameters.ChannelNumPar = Npixel;
gParas_Inside.VariableParameters.ChannelParSpace = prmflow.recon.delta_d;
% gParas_Inside.ReconParameters.reconKernel = FBPprm.reconKernel;
gParas_Inside.ReconParameters.reconKernel = reconKernel;

% prepare the convolution
subConv = Conv('Conv', 'GlobalParameter', gParas_Inside, 'DebugSave', 0, 'ProcessType', 'Matlab');
subConv.Init(gParas_Inside);

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
    BpStruct.nSliceNumber        = 1;                
    BpStruct.fSliceSpace         = prmflow.system.detector.hz_ISO;
    % BpStruct.nCouchDirection     = gParas.AcquisitionParameter.tableDirection; % Add for 3DBP
    % BpStruct.fSourseToIsocenter  = gParas.GeometryParameter.sourceToIso; % Add for 3DBP
    BpStruct.fTiltAngle          = prmflow.protocol.gantrytilt;
    BpStruct.fZDFS               = 0; % zDFS 情况下光源抖动偶数View位置相对奇数View位置的距离
    BpStruct.fMaxFOV             = 500;
    BpStruct.fReconFOV           = reconFOV;
    BpStruct.nXPixels            = imagesize;
    BpStruct.nYPixels            = imagesize;
    BpStruct.fXReconCenter       = 0;
    BpStruct.fYReconCenter       = 0;
    % BpStruct.fImageThickness     = gParas.ReconParameters.ImageThickness; % Add for 3DBP, NO Using
%     BpStruct.fImageIncrement     = prmflow.recon.imageincrement; % Add for 3DBP
    BpStruct.nImageNumber        = 1;

    fViewWeight = ones(BpStruct.nTotalViewNumber, 1, 'single')./4;
    ImgNum_pZ = BpStruct.nSliceNumber;
    
    % QDO fix
    if isQDO
        struct_raw_Inside.data = cat(3, struct_raw_Inside.data, zeros(size(struct_raw_Inside.data)));
        struct_raw_Inside.raw_size(3) = Nviewprot*2;
        BpStruct.nViewPerRevolution  = Nviewprot*2;
        BpStruct.nTotalViewNumber    = Nviewprot*2;
        fViewWeight = fViewWeight.*2;
    end
    
    func = Engine.GetFunction('Axial2DBP_CU');
    for islice = 1:Nslice
        img_index = Nslice*(ishot-1) + islice;
        BpStruct.fXReconCenter       = -prmflow.recon.imagecenter(img_index, 2);
        BpStruct.fYReconCenter       = prmflow.recon.imagecenter(img_index, 1);
        % I know x is y, y is x. Y? yeX
        CorrImage = func(squeeze(struct_raw_Inside.data(:,islice,:)), fViewWeight, BpStruct, 1);
        dataflow.image(:,:,img_index) = CorrImage';
    end
end

% reorder
reorderflag = prmflow.protocol.couchdirection < 0;
dataflow.image = imagereorder(dataflow.image, Nslice, reorderflag);

%% done
% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end