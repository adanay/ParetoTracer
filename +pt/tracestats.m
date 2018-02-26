% Copyright (c) 2018 Adanay Martín & Oliver Schütze.
% This file is subject to the terms and conditions defined in
% the file 'LICENSE.txt', which is part of this source code package.

function [stats] = tracestats(stats, optstats)
% If used without arguments, initializes the Pareto Tracer method 
% statistics structure.
% Otherwise, updates the PT statistics values with those obtained by the 
% corrector phase. 

if nargin == 0
  stats = struct(...
    'PCIts', 0,...
    'OptIts', 0,...
    'OptDirIts', 0,...
    'OptLsIts', 0,...
    'OptCount',0,...
    'LastPCIt', 0,...
    'LastOptIts', 0,...
    'LastOptDirIts', 0,...
    'LastOptLsIts', 0,...
    'LastOptCount', 0,...
    'fCount', 0,...
    'JCount', 0,...
    'HCount', 0,...
    'JvCount', 0, 'wJCount', 0, 'wJvCount', 0, 'vHCount', 0, 'HwCount', 0, 'HwvCount', 0,...
    'aCount', 0,...
    'aeqCount', 0,...
    'cCount', 0,...
    'JcCount', 0,...
    'HcCount', 0,...
    'JcvCount', 0, 'wJcCount', 0, 'wJcvCount', 0, 'vHcCount', 0, 'HcwCount', 0, 'HcwvCount', 0,...
    'ceqCount', 0,...
    'JceqCount', 0,...
    'HceqCount', 0,...
    'JceqvCount', 0, 'wJceqCount', 0, 'wJceqvCount', 0, 'vHceqCount', 0, 'HceqwCount', 0, 'HceqwvCount', 0,...
    'Count', 0);
else
  [stats] = updstats(stats, optstats);
end
end

function [stats] = updstats(stats, optstats)
% Update the Pareto Tracer statistics with those obtained by the corrector
% phase.

stats.OptIts = stats.OptIts + optstats.OptIts;
stats.OptDirIts = stats.OptDirIts + optstats.OptDirIts;
stats.OptLsIts = stats.OptLsIts + optstats.OptLsIts;
stats.OptCount = stats.OptCount + 1;

if stats.LastPCIt == stats.PCIts
  stats.LastOptIts = stats.LastOptIts + optstats.OptIts;
  stats.LastOptDirIts = stats.LastOptDirIts + optstats.OptDirIts;
  stats.LastOptLsIts = stats.LastOptLsIts + optstats.OptLsIts;
  stats.LastOptCount = stats.LastOptCount + 1;
else
  stats.LastPCIt = stats.PCIts;
  stats.LastOptIts = optstats.OptIts;
  stats.LastOptDirIts = optstats.OptDirIts;
  stats.LastOptLsIts = optstats.OptLsIts;
  stats.LastOptCount = 1;
end

stats.fCount = stats.fCount + optstats.fCount;
stats.JCount = stats.JCount + optstats.JCount;
stats.HCount = stats.HCount + optstats.HCount;
stats.JvCount = stats.JvCount + optstats.JvCount; 
stats.wJCount = stats.wJCount + optstats.wJCount;  
stats.wJvCount = stats.wJvCount + optstats.wJvCount;  
stats.vHCount = stats.vHCount + optstats.vHCount;  
stats.HwCount = stats.HwCount + optstats.HwCount;  
stats.HwvCount = stats.HwvCount + optstats.HwvCount; 
stats.aCount = stats.aCount + optstats.aCount; 
stats.aeqCount = stats.aeqCount + optstats.aeqCount;
stats.cCount = stats.cCount + optstats.cCount; 
stats.JcCount = stats.JcCount + optstats.JcCount; 
stats.HcCount = stats.HcCount + optstats.HcCount;
stats.JcvCount = stats.JcvCount + optstats.JcvCount;  
stats.wJcCount = stats.wJcCount + optstats.wJcCount;  
stats.wJcvCount = stats.wJcvCount + optstats.wJcvCount;  
stats.vHcCount = stats.vHcCount + optstats.vHcCount;  
stats.HcwCount = stats.HcwCount + optstats.HcwCount;  
stats.HcwvCount = stats.HcwvCount + optstats.HcwvCount; 
stats.ceqCount = stats.ceqCount + optstats.ceqCount; 
stats.JceqCount = stats.JceqCount + optstats.JceqCount; 
stats.HceqCount = stats.HceqCount + optstats.HceqCount; 
stats.JceqvCount = stats.JceqvCount + optstats.JceqvCount;  
stats.wJceqCount = stats.wJceqCount + optstats.wJceqCount; 
stats.wJceqvCount = stats.wJceqvCount + optstats.wJceqvCount;  
stats.vHceqCount = stats.vHceqCount + optstats.vHceqCount;  
stats.HceqwCount = stats.HceqwCount + optstats.HceqwCount;  
stats.HceqwvCount = stats.HceqwvCount + optstats.HceqwvCount; 
end

