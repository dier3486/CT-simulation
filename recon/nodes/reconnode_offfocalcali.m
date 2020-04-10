function [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status)

% generate the offfocal as a node 
% 20200410


offfocalpara = myxml2struct('E:\PANGU.DAT\CT\SINO\Config\bhcali\offfocalparameter.xml');
% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nps = Npixel*Nslice;
Nview = prmflow.recon.Nview;

% parameters to use
caliprm = prmflow.pipe.(status.nodename);

% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end

xml = dataflow.pdxml;
structxml = Global_Parameter(xml);
% structxml.RawData.CortblPath = CorFileDir;
CorFindStruct = GetCorFindStruct(structxml);

KVp = CorFindStruct.kV;
Bowtie = structxml.AcquisitionParameter.bowTie;
collimator = CorFindStruct.Collimation;
singlecollimtorwidth = CorFindStruct.SliceThickness;

Num = size(offfocalpara.configure.offfocalpara,2);
    
for jj = 1:Num
    kVpara = offfocalpara.configure.offfocalpara{1, jj}.KVp;
    Bowtiepara = offfocalpara.configure.offfocalpara{1, jj}.Bowtietype;
    if kVpara == KVp && strcmp(Bowtiepara,Bowtie)
        Nmcollimation = size(offfocalpara.configure.offfocalpara{1, jj}.collimation,2);
        for kk = 1:Nmcollimation
            Collimationpara = offfocalpara.configure.offfocalpara{1, jj}.collimation{1, kk}.collimationwidth;
            Acqcollimation = [num2str(collimator) 'x' num2str(singlecollimtorwidth)];
            if strcmp(Collimationpara,Acqcollimation)
                offintensity = offfocalpara.configure.offfocalpara{1, jj}.collimation{1, kk}.offintensity;
                offwidth = offfocalpara.configure.offfocalpara{1, jj}.collimation{1, kk}.offwidth;
                offedge = offfocalpara.configure.offfocalpara{1, jj}.collimation{1, kk}.offedge;
                ratescale = offfocalpara.configure.offfocalpara{1, jj}.collimation{1, kk}.ratescale;
                break;
            end
            
        end
        break;
    end
end


% paramters for corr
offfocalcorr = caliprmforcorr(prmflow, corrversion);
offfocalcorr.offintensity = offintensity;
offfocalcorr.offwidth = offwidth;
offfocalcorr.offedge = offedge;
offfocalcorr.ratescale = ratescale;
offfocalcorr.order = 1;
offfocalcorr.ratesize = dataflow.beamhardencorr.ratesize;
offfocalcorr.main = dataflow.beamhardencorr.airrate;

% to return 
dataflow.offfocalcorr =  offfocalcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end