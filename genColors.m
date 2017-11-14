function colors = genColors(n, style, range)
% Generate sequential colors for plotting lines in order.
%
% Inputs:
%   n [ scalar positive integer ]
%       Numer of color combos to generate
%   style [ {'PuBu'} | 'OrRd' ]
%       Style/colormap to use
%   start [ 1 x 2 double vector {[0,1]} ]
%       Optional range of colors from colormap. Set the first value to > 0 if
%       the starting colors are too light.
%
% Outputs:
%   colors [ n x 3 double matrix ]
%       Ouput colors interpolated from style/colormap in RGB format
%
% Source: http://colorbrewer2.org/
%
% TODO: implement other styles

if nargin < 3
    range = [];
    if nargin < 2
        style = [];
    end
end

if isempty(style)
    style = 'PuBu';
end
if isempty(range)
    range = [0,1];
end

assert(all(size(range) == [1,2]), 'genColors:InvalidRangeDims', 'range must be a 1x2 double vector')
assert(isnumeric(range) && all(range >= 0) && all(range <= 1) && range(1) <= range(2), 'genColors:InvalidRange', 'range must be numbers between 0 and 1, where 2nd val is greater than first')

%% Select colormap
switch style

    
    % 9-class PuBu/darkest 8
    %   Sequential coloring
    case 'PuBu'
        cmap = [236,231,242;
            208,209,230;
            166,189,219;
            116,169,207;
            54,144,192;
            5,112,176;
            4,90,141;
            2,56,88];
        
    % 9-class PuBu
    % cmap = [255,247,251; % first few lines are too light
    %     236,231,242;
    %     208,209,230;
    %     166,189,219;
    %     116,169,207;
    %     54,144,192;
    %     5,112,176;
    %     4,90,141;
    %     2,56,88]
    
    % 9-class OrRd/darkest 8
    %   Sequential coloring
    case 'OrRd'
        cmap = [254,232,200;
            253,212,158;
            253,187,132;
            252,141,89;
            239,101,72;
            215,48,31;
            179,0,0;
            127,0,0];
        
%         cmap = [255,247,236; % first few lines are too light
%             254,232,200;
%             253,212,158;
%             253,187,132;
%             252,141,89;
%             239,101,72;
%             215,48,31;
%             179,0,0;
%             127,0,0];
        
    otherwise
        error('Style %s not recognized', style)
end

%% Normalize colormap
% cmapMax = max(max(cmap));
cmapMax = 255; % max value of uint32
cmap = cmap ./ cmapMax;
x = linspace(0,1,size(cmap,1))';

%% Get colors
colors = zeros(n,3);
xq = linspace(range(1), range(2), n)';
for i = 1:3
    colors(:,i) = interp1(x, cmap(:,i), xq);
end