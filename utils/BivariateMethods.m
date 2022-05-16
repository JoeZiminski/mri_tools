
classdef BivariateMethods
    % Methods for correlation and regression
    
    % Requires most recent github version of Robust Correlation Toolbox to
    % calculate Pearson and Skipped Correlation. 
    % https://github.com/CPernet/Robust-Correlations
    
    % jjz33@cam.ac.uk
    % ---------------------------------------------------------------------------------------------------------------
    
    properties
        
        plotter = Plotter;    
        
        
    end
    
    
    methods (Static)
        
        function [r, p, betas, weights] = robust_regression(x, y, x_label, y_label)
            %
            % Calculate robust regression. For the r value take the Rsquared from the model
            % output, sqrt and multiply by the sign of the t statistic for the slope.
            %
            % Take the p value from the coefficient (not the model!)
            %
            % See self.correlate_package() for inputs. Outputs parameters + prints results.
            %
            % NOTE: Only good for 1 comparison as takes the p value for the x1 coefficient)
            % --------------------------------------------------------------------------
            
            if size(x, 2) > 1
                warning('p value returned is only for first predictor. This function is for use',  ...
                        'with one x / y pair only');
            end
            
            fitlm_mdl = fitlm(x, y, 'RobustOpts', 'on');
            r = sqrt(abs(fitlm_mdl.Rsquared.Ordinary)) .* sign(fitlm_mdl.Coefficients{2, 'tStat'});
            p = fitlm_mdl.Coefficients{'x1', 'pValue'};
            
            fprintf('robust_correlation: r: %.2f, p: %.3f, for %s vs %s\n', r, p, x_label, y_label);

            betas = table(fitlm_mdl.Coefficients{1, 1}, fitlm_mdl.Coefficients{2, 1}, 'VariableNames', {'intercept', 'slope'});

            weights = fitlm_mdl.Robust.Weights;
            
        end
        
       
        function [r, p, hboot, CI] = pearson(x, y, x_label, y_label)
            % Perform Pearson correlation, print result to command window
            % and return r / sig.
            % Robust Correlation Toolbox must be on path.
            % ------------------------------------------------------------         
            [r, ~, p, hboot, CI] = Pearson(x, y, 0);
            fprintf('Pearson: r: %.3f, : %.3f, hboot %1.0f CI: [%.2f, %.2f] for %s vs %s\n', r, p, hboot, CI, x_label, y_label);
        
        end
        
        
        function [r, sig, outids] = skipped_corr(x, y, x_label, y_label)
            % Skipped correlate as implimented by the Robust Correlation
            % Toolbox (github version). Prints result and returns in
            % structure.
            % -----------------------------------------------------------
            [r, ~, h, outid] = run_skipped_correlation(x, y, 0);
            fprintf('skipped correlation: r: %.2f, h: %1.0f, outliers: %d for %s vs %s\n', r.Pearson, h.Pearson, length(outid{:}), x_label, y_label);
 
            r = r.Pearson;
            if h.Pearson == 1  
                sig = 'sig';
            else
                sig = 'n.s.';       
            end
            if isempty(outid{:}) 
                outids = 0;
            else
                outids = length(outid{:});
            end

        end
        
        
        function [resid] = regress_data(x, y, robust)
            % Regress data X from Y, returning residuals. 
            %
            % if robust = 'robust', robust regression is used.
            % ----------------------------------------------------------
            [~, p] = corr(x, y);
            if p < 0.05
                warning(' Regressing one variable from another - the two variables are correlated')
            end

            if strcmp(robust, 'robust')
                [~, stats] = robustfit(x, y);
                resid = stats.resid;
            elseif strcmp(robust, 'not_robust')
                [~, ~, resid] = regress(y, [ones(length(x), 1), x]);
            else
                error('input must be ''robust'' or ''not_robust''');
                
            end
            
        end
        
        function check_univariate_outlier(data_table) 

            isoutlier_bool = isoutlier(data_table{:, 1}, 'grubbs');

            if any(isoutlier_bool)
                outlier_id = data_table.Properties.RowNames(isoutlier_bool);

                for id = outlier_id
                    fprintf("outlier detected: %s\n", char(id));
                end

            else

                fprintf("no outliers detected\n");

            end

        end
        
    end
    
end
