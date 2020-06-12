% pin cali main, test
pinxml = 'F:\data-Dier.Z\PG\bay3_0601\Pinrecon.xml';
pindata = 'F:\data-Dier.Z\PG\bay3_0601\pin\pin_1';
[~, dataflow_all, prmflow_all] = CRISrecon(pinxml, pindata);

Nview = prmflow_all.protocol.viewnumber;
Nviewprot = prmflow_all.protocol.viewperrot;

viewstep = Nviewprot/4;
Nstep = 20;
% x0 = [200 0.2 0 0 1 0 1 0 0];

prmflow = prmflow_all;
prmflow.recon.Nview = Nviewprot;

figid = figure;
fanfix = cell(1, Nstep);
for istep = 1:Nstep
    startview = (istep-1)*viewstep + 1;
    endview = startview + Nviewprot - 1;
    if endview > Nview
        break;
    end
    dataflow = dataflow_all;
    dataflow.rawdata = dataflow.rawdata(:, startview:endview);
    for ifield = fieldnames(dataflow.rawhead)'
        dataflow.rawhead.(ifield{1}) = dataflow.rawhead.(ifield{1})(:, startview:endview);
    end
    fanfix{istep} = pindetcali_test1(dataflow, prmflow);
    figure(figid);
    plot(mean(fanfix{istep}, 2));
    drawnow;
end


