function [] = save_figure(fighandle,folder_name)
%SAVE_PATCHFIG to save patch figures with high resolution
 saveas(fighandle,strcat(folder_name,'.eps'),'epsc')
 saveas(fighandle,strcat(folder_name,'.png'),'png')
 saveas(fighandle,strcat(folder_name,'.svg'),'svg');
 saveas(fighandle,strcat(folder_name,'.fig'),'fig');
end

