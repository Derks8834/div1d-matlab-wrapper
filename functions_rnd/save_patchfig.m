function [] = save_patchfig(fighandle,folder_name,varargin)

if size(varargin) < 1
    resolution = 1200;
else
    resolution = varargin;
end
%SAVE_PATCHFIG to save patch figures with high resolution
 exportgraphics(fighandle,strcat(folder_name,'.eps'),'Resolution',resolution)
 exportgraphics(fighandle,strcat(folder_name,'.png'),'Resolution',resolution)
 saveas(fighandle,strcat(folder_name,'.svg'),'svg');
 saveas(fighandle,strcat(folder_name,'.fig'),'fig');
end

