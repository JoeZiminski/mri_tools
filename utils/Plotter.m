classdef Plotter
    
    methods 
        
        function [self] = Plotter()
           % Constructer
        end


        function make_heatmap(~, data, labels_x, labels_y, title_, caxis_, varargin)
            % Make a heatmap for correlation matricies
            %
            % INPUT:  data - 2D matrix for plotting
            %         labels - cell array of axis labels e.g. {'V1', 'LOC'}
            %         caxis - scale of the colour 'axis' e.g. [0, 1]
            %         varargin - colormap to use (default 'jet')
            %
            % -------------------------------------------------------------
            
            % remove underscore from labels as useless matlab cannot set interpreter on the heatmap
            for i = 1:length(labels_x)
                labels_x{i} = strrep(labels_x{i}, '_', ' ');
            end
            for i = 1:length(labels_y)
                labels_y{i} = strrep(labels_y{i}, '_', ' ');
            end
            
            % set default color 
            if ~isempty(varargin)
                colour = varargin{1};
            else
                colour = 'jet';
            end
            
            % plot 
            figure('Position', [500, 300, 800, 600], 'Color', 'white');
            hm = heatmap(data, 'XDisplayLabels', labels_x, 'YDisplayLabels', labels_y);  % 'CellLabelColor', 'none'
            hm.CellLabelFormat = '%0.2g';
            hm.FontSize = 15;

            % put x axis on top
            axp = struct(hm);  % can ignore warning
            axp.Axes.XAxisLocation = 'top';
            
            colormap(colour);
            caxis(caxis_);
            title_ = strrep(title_, '_', '\_');  % heatmap only supports tex? formatting
            title(title_)% moving title from top is an absolute nightmare - why is matlab plotting API so bad?

        end

        function make_weighted_scatterplot(~, x, y, x_label, y_label,out) 
           %
           % Make a scatterplot
           %
           % TODO: weighting option removed, rename 
           % --------------------------------------------------------------------------
           f = figure('Color', 'white', 'Position', [500, 300, 750, 600]);  
           num_points = length(y);
           line_width = 1.5; 
           font_size = 14;

           green = [0/255, 153/255, 0/255];
           blue = [0/255, 104/255, 239/255];
           red = [226/255, 19/255, 3/255];
           color = blue;     

           ax = axes;
           hold on;
           for i = 1:length(x)
                s = scatter(x(i), y(i), 240, color, 'filled');
           end 

           xlabel(x_label, 'Interpreter', 'None');
           ylabel(y_label, 'Interpreter', 'None');
           x_pad = range(x) * 0.1;
           if isfield(out, 'betas')
               intercept =  out.('betas').intercept;
               slope =  out.('betas').slope;
               x_fit = linspace(min(x) - x_pad, max(x) + x_pad, 10000);
               fit = intercept + x_fit * slope;
               hold on; plot(x_fit, fit, 'LineWidth', 6, 'Color', color); 
           end

           ax.LineWidth = line_width;
           ax.FontSize = font_size;

        end

        function [stderr] = calculate_stderr(~, data)
        % Input must be N x p and will return 1 x p array of
        % standard errors
        % ---------------------------------------------
            stderr = std(data)/sqrt(size(data,  1));
        end
    
        function make_bar_graph(self, data, bar_color, errorbar_color, xticklabels_, plot_title, ylim_, norm_data_flag, varargin)
            % Make a 2 - bar graph where cols 2-3 are normalised
            % to column 1.
            %
            % varargin - optionally give position for figure
            % --------------------------------------------------

            line_width = 1.5;
            
            if istable(data)
                data = table2array(data);
            end
            
            if bar_color == false
                bar_color = [0 0.4470 0.7410];
            end
                 
            if norm_data_flag
                % express as a percentage of first session
                data_norm = (data - data(:, 1)) ./ data(:, 1); 
                data_norm = data_norm * 100;
                data = data_norm(:, 2:3);  
                
            end 
                       
            if ~isempty(varargin) 
                position = varargin{1};
                f = figure('Color', 'white', 'Position', position);
            else
                f = figure('Color', 'white');
            end
            
            bar_ = bar(mean(data), 0.7, 'FaceColor', bar_color, 'LineWidth', 1); hold on;
            bar_.EdgeColor = 'none';

            e = errorbar(mean(data), self.calculate_stderr(data),  'LineWidth',line_width, 'Color', errorbar_color);  % 'Color', errorbar_color,
            e.LineStyle = 'none';
            
            set(gca, 'box', 'off', 'LineWidth', line_width, 'FontSize', 14);
            xticklabels(xticklabels_);
            title(plot_title, 'Interpreter', 'None');
            
            if ylim_
                ylim(ylim_);
            end
            
            set(gca, 'TickLength',  [0, 0]);

        end
               
    end 
end
                  
