function label = makeLabelFromCereConfig(config,pair)
vals = [config.Config];
ppair = pair+1;
ma = str2num(vals{3})/1000;
pw = str2num(vals{5});
f = str2num(vals{7});
label = sprintf("Pair %d-%d, %d ma, %d %ss, %d Hz", ...
    pair,ppair,ma,pw,char(181),f);
end