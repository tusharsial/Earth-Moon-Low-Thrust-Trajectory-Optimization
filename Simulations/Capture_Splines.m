function [vr_spline, vtheta_spline, tf_spline] = Capture_Splines(rf_vec, m_llo_vec, tf_vec, vr_vec, vtheta_vec)
    % CREATE_VELOCITY_SPLINES Creates 2D cubic spline interpolants
    % 
    % Input organization:
    %   tf_vec, vr_vec, vtheta_vec are 2D arrays with:
    %   - rows corresponding to m_llo_vec (spacecraft mass)
    %   - columns corresponding to rf_vec (moon radial distance)
    %
    % Outputs:
    %   Three griddedInterpolant objects in NDGRID format

    % Verify input dimensions
    [n_rows, n_cols] = size(tf_vec);
    if ~isequal(size(vr_vec), [n_rows, n_cols]) || ~isequal(size(vtheta_vec), [n_rows, n_cols])
        error('All input 2D arrays must have the same dimensions');
    end
    
    if length(m_llo_vec) ~= n_rows || length(rf_vec) ~= n_cols
        error('Dimension mismatch: length(m_llo_vec) must match rows, length(rf_vec) must match columns');
    end

    % Convert to NDGRID format that griddedInterpolant expects
    % Note the order of inputs to ndgrid matches your data organization
    [rf_grid_nd, m_llo_grid_nd] = ndgrid(rf_vec, m_llo_vec);
    
    % Create the interpolants - transposing the data arrays is CRUCIAL
    vr_spline = griddedInterpolant(rf_grid_nd, m_llo_grid_nd, vr_vec', 'spline');
    vtheta_spline = griddedInterpolant(rf_grid_nd, m_llo_grid_nd, vtheta_vec', 'spline');
    tf_spline = griddedInterpolant(rf_grid_nd, m_llo_grid_nd, tf_vec', 'spline');
end