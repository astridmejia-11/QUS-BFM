
function save_all_figures_to_directory(dir_name,title, varargin)
% optional argument is other format to save the image

title = string(title);

figlist=findobj('type','figure');
number=get(figlist,'Number');
for i=1:numel(figlist)

    fileName = fullfile(dir_name,title+number(i)+".png");
    exportgraphics(figlist(i), fileName, 'Resolution','300')    
    if nargin == 3
        fileName = fullfile(dir_name,title+number(i)+"."+varargin{1});
        saveas(figlist(i), fileName);
    end
    pause(0.1)
end

end