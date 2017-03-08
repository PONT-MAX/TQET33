
% Qins.m : Insert a new element into the SNV queue
%
% snvQout = Qins(snvQin,snv,topnode)
%
% snvQin:   SNV queue with SNVs in column 1, top nodes in column 2
% topnode:  topnode of the element to be inserted
% snv:      SNV of the element to be inserted
%
% snvQout:  updated SNV queue
%
% RGB INI September 1998

function snvQout = Qins(snvQin,snv,topnode);
snvQout = [snvQin; [topnode,snv]];
return



